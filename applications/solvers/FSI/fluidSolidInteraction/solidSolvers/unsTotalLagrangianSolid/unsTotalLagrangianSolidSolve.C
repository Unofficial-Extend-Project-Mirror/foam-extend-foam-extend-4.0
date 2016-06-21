/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "unsTotalLagrangianSolid.H"

#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGradf.H"

#include "tractionDisplacementFvPatchVectorField.H"
#include "skewCorrectionVectors.H"
#include "multiMaterial.H"
#include "twoDPointCorrector.H"

#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidSolvers
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volVectorField> unsTotalLagrangianSolid::residual
(
    const volVectorField& source
)
{
    Switch nonLinear(solidProperties().lookup("nonLinear"));

    dimensionedScalar K("K", dimless/dimTime, 0);
    if (solidProperties().found("K"))
    {
        K = dimensionedScalar(solidProperties().lookup("K"));
    }

    fvVectorMatrix DEqn
    (
        rho_*fvm::d2dt2(D_)
      - fvm::laplacian(2*muf_ + lambdaf_, D_, "laplacian(DD,D)")
     == fvc::div
        (
            mesh().Sf()
          & (
              - (muf_ + lambdaf_)*gradDf_
              + muf_*gradDf_.T() + lambdaf_*(I*tr(gradDf_))
            )
        )
    );

    // Add damping
    if (K.value() > SMALL)
    {
        DEqn += K*rho_*fvm::ddt(D_);
    }

    if (nonLinear)
    {
        surfaceSymmTensorField Ef =
            symm(gradDf_) + 0.5*symm(gradDf_ & gradDf_.T());

        surfaceSymmTensorField sigmaf = 2*muf_*Ef  + I*(lambdaf_*tr(Ef));

        DEqn -=
            fvc::div
            (
                muf_*(mesh().Sf() & (gradDf_ & gradDf_.T()))
              + 0.5*lambdaf_*tr(gradDf_ & gradDf_.T())*mesh().Sf()
            )
          + fvc::div(mesh().Sf() & (sigmaf & gradDf_));
    }

    if (interface().valid())
    {
        interface()->correct(DEqn);
    }

    DEqn += source;

    tmp<volVectorField> tResidual
    (
        new volVectorField
        (
            IOobject
            (
                "residual(" + D_.name() + ")",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedVector("0", dimForce, vector::zero)
        )
    );

    tResidual().internalField() = DEqn.residual();


//     // Boundary residual

//     D_.storePrevIter();
//     D_.boundaryField().evaluate();

//     forAll(tResidual().boundaryField(), patchI)
//     {
//         tResidual().boundaryField()[patchI] =
//             D_.prevIter().boundaryField()[patchI]
//           - D_.boundaryField()[patchI];
//     }

//     D_ = D_.prevIter();

    return tResidual;
}


void unsTotalLagrangianSolid::initialise
(
    const volVectorField& sol,
    bool consistentBoundaryField
)
{
    Info << "Initialise solid solver solution" << endl;

    D_ = sol;

    if (interface().valid())
    {
        interface()->updateDisplacement(pointD_);
        interface()->updateDisplacementGradient(gradD_, gradDf_);
    }
    else
    {
        volToPoint_.interpolate(D_, pointD_);
        gradD_ = fvc::grad(D_, pointD_);
        gradDf_ = fvc::fGrad(D_, pointD_);
    }

    D_.correctBoundaryConditions();

    // Correct boundary conditions
    if (consistentBoundaryField)
    {
        PtrList<vectorField> prevD(D_.boundaryField().size());
        forAll(D_.boundaryField(), patchI)
        {
            prevD.set
            (
                patchI,
                new vectorField
                (
                    D_.boundaryField()[patchI].size(),
                    vector::zero
                )
            );
        }

        scalar residual = GREAT;
        label iCorr = 0;
        do
        {
            // Recalculate boundary points displacement
            forAll (pointD_.boundaryField(), patchI)
            {
                vectorField patchPointU =
                    volToPoint_.interpolate(mesh().boundaryMesh()[patchI], D_);

                pointD_.boundaryField()[patchI].setInInternalField
                (
                    pointD_.internalField(),
                    patchPointU
                );
            }
            pointD_.correctBoundaryConditions();


            // Recalculate boundary face gradient
            forAll(mesh().boundary(), patchI)
            {
                if (mesh().boundary()[patchI].size())
                {
                    vectorField pPointD =
                        pointD_.boundaryField()[patchI].patchInternalField();

                    tensorField pGradDf =
                        fvc::fGrad(mesh().boundaryMesh()[patchI], pPointD);

                    vectorField n = mesh().boundary()[patchI].nf();

                    gradDf_.boundaryField()[patchI] =
                        pGradDf + n*D_.boundaryField()[patchI].snGrad();
                }
            }

            // Recalculate boundary cell gradient
            // (may be not)

            // Save previous boundary displacement
            forAll(D_.boundaryField(), patchI)
            {
                if (D_.boundaryField().size())
                {
                    prevD[patchI] = D_.boundaryField()[patchI];
                }
            }

            // Correct boundary conditions
            D_.correctBoundaryConditions();

            // Calculate relative residual
            scalar maxD = 0;
            forAll(D_.boundaryField(), patchI)
            {
                if (D_.boundaryField()[patchI].size())
                {
                    scalar maxPatchD = max(mag(D_.boundaryField()[patchI]));

                    if (maxPatchD > maxD)
                    {
                        maxD = maxPatchD;
                    }
                }
            }
            residual = 0;
            forAll(D_.boundaryField(), patchI)
            {
                if (D_.boundaryField()[patchI].size())
                {
                    scalar maxPatchResidual =
                        max
                        (
                            mag(D_.boundaryField()[patchI] - prevD[patchI])
                           /(maxD + SMALL)
                        );

                    if (maxPatchResidual > residual)
                    {
                        residual = maxPatchResidual;
                    }
                }
            }
        }
        while ( (residual > 1e-10) && (++iCorr<=50) );
//         while (++iCorr<=50);

        gradDf_ = fvc::fGrad(D_, pointD_);

        Info << "iCorr: " << iCorr
            << ", max rel residual: " << residual << endl;
    }
}


scalar unsTotalLagrangianSolid::smooth
(
    const volVectorField& source,
    label nCorrectors
)
{
    Info << "Smoothing solid solver solution" << endl;

    label nCorr = nCorrectors;
    if (nCorrectors == 0)
    {
        nCorr = readInt(solidProperties().lookup("nCorrectors"));
    }

    Switch nonLinear(solidProperties().lookup("nonLinear"));

    Switch debug(solidProperties().lookup("debug"));

    dimensionedScalar K("K", dimless/dimTime, 0);
    if (solidProperties().found("K"))
    {
        K = dimensionedScalar(solidProperties().lookup("K"));
    }

    int iCorr = 0;
    scalar initialResidual = 0;
    lduSolverPerformance solverPerf;
    scalar residual = GREAT;

    lduMatrix::debug = debug;

    do
    {
        if (lduMatrix::debug)
        {
            Info<< "Time: " << runTime().timeName()
                << ", outer iteration: " << iCorr << endl;
        }

        D_.storePrevIter();

        fvVectorMatrix DEqn
        (
            rho_*fvm::d2dt2(D_)
          - fvm::laplacian(2*muf_ + lambdaf_, D_, "laplacian(DD,D)")
         == fvc::div
            (
                mesh().Sf()
              & (
                  - (muf_ + lambdaf_)*gradDf_
                  + muf_*gradDf_.T() + lambdaf_*(I*tr(gradDf_))
                )
            )
        );

        // Add damping
        if (K.value() > SMALL)
        {
            DEqn += K*rho_*fvm::ddt(D_);
        }

        if (nonLinear)
        {
            surfaceSymmTensorField Ef =
                symm(gradDf_) + 0.5*symm(gradDf_ & gradDf_.T());

            surfaceSymmTensorField sigmaf = 2*muf_*Ef  + I*(lambdaf_*tr(Ef));

            DEqn -=
                fvc::div
                (
                    muf_*(mesh().Sf() & (gradDf_ & gradDf_.T()))
                  + 0.5*lambdaf_*tr(gradDf_ & gradDf_.T())*mesh().Sf()
                )
              + fvc::div(mesh().Sf() & (sigmaf & gradDf_));
        }

        if (interface().valid())
        {
            interface()->correct(DEqn);
        }

        DEqn += source;

        solverPerf = DEqn.solve();

//         forAll(D_.boundaryField(), patchI)
//         {
//             D_.boundaryField()[patchI] += source.boundaryField()[patchI];
//         }

        if (iCorr == 0)
        {
            initialResidual = solverPerf.initialResidual();
        }

        D_.relax();

        if (interface().valid())
        {
            interface()->updateDisplacement(pointD_);
            interface()->updateDisplacementGradient(gradD_, gradDf_);
        }
        else
        {
            volToPoint_.interpolate(D_, pointD_);
            gradD_ = fvc::grad(D_, pointD_);
            gradDf_ = fvc::fGrad(D_, pointD_);
        }


        // Calculate momentum residual
        {
            scalar maxDD =
                gMax
                (
                    mag(D_.internalField() - D_.oldTime().internalField())
                );

            residual =
                gMax
                (
                    mag(D_.internalField() - D_.prevIter().internalField())
                   /(maxDD + SMALL)
                );

            if (lduMatrix::debug)
            {
                Info << "Relative residual = " << residual << endl;
            }
        }
    }
    while(++iCorr < nCorr);

    U_ = fvc::ddt(D_);

    // Calculate second Piola-Kirchhoff stress
    {
        volSymmTensorField E = symm(gradD_);

        if(nonLinear)
        {
            E += 0.5*symm(gradD_ & gradD_.T());
        }

        sigma_ = 2*mu_*E + I*(lambda_*tr(E));
    }

    Info << solverPerf.solverName() << ": Solving for " << D_.name()
        << ", Initial residula = " << initialResidual
        << ", Final residual = " << solverPerf.initialResidual()
        << ", No outer iterations = " << iCorr
        << ", Relative momentum residual = " << residual << endl;

    lduMatrix::debug = 1;

    return residual;
}


scalar unsTotalLagrangianSolid::smooth
(
    const volVectorField& source,
    const volVectorField& refSolution,
    label nCorrectors
)
{
    Info << "Smoothing solid solver solution" << endl;

    label nCorr = nCorrectors;
    if (nCorrectors == 0)
    {
        nCorr = readInt(solidProperties().lookup("nCorrectors"));
    }

    Switch nonLinear(solidProperties().lookup("nonLinear"));

    Switch debug(solidProperties().lookup("debug"));

    dimensionedScalar K("K", dimless/dimTime, 0);
    if (solidProperties().found("K"))
    {
        K = dimensionedScalar(solidProperties().lookup("K"));
    }

    int iCorr = 0;
    scalar initialResidual = 0;
    lduSolverPerformance solverPerf;
    scalar residual = GREAT;

    lduMatrix::debug = debug;

    do
    {
        if (lduMatrix::debug)
        {
            Info<< "Time: " << runTime().timeName()
                << ", outer iteration: " << iCorr << endl;
        }

        D_.storePrevIter();

        fvVectorMatrix DEqn
        (
            rho_*fvm::d2dt2(D_)
          - fvm::laplacian(2*muf_ + lambdaf_, D_, "laplacian(DD,D)")
         == fvc::div
            (
                mesh().Sf()
              & (
                  - (muf_ + lambdaf_)*gradDf_
                  + muf_*gradDf_.T() + lambdaf_*(I*tr(gradDf_))
                )
            )
        );

        // Add damping
        if (K.value() > SMALL)
        {
            DEqn += K*rho_*fvm::ddt(D_);
        }

        if (nonLinear)
        {
            surfaceSymmTensorField Ef =
                symm(gradDf_) + 0.5*symm(gradDf_ & gradDf_.T());

            surfaceSymmTensorField sigmaf = 2*muf_*Ef  + I*(lambdaf_*tr(Ef));

            DEqn -=
                fvc::div
                (
                    muf_*(mesh().Sf() & (gradDf_ & gradDf_.T()))
                  + 0.5*lambdaf_*tr(gradDf_ & gradDf_.T())*mesh().Sf()
                )
              + fvc::div(mesh().Sf() & (sigmaf & gradDf_));
        }

        if (interface().valid())
        {
            interface()->correct(DEqn);
        }

        DEqn += source;

        solverPerf = DEqn.solve();

        if (iCorr == 0)
        {
            initialResidual = solverPerf.initialResidual();
        }

        D_.relax();

        if (interface().valid())
        {
            interface()->updateDisplacement(pointD_);
            interface()->updateDisplacementGradient(gradD_, gradDf_);
        }
        else
        {
            volToPoint_.interpolate(D_, pointD_);
            gradD_ = fvc::grad(D_, pointD_);
            gradDf_ = fvc::fGrad(D_, pointD_);
        }


        // Calculate momentum residual
        {
            scalar maxDD =
                gMax
                (
                    mag(D_.internalField() - D_.oldTime().internalField())
                );

            residual =
                gMax
                (
                    mag(D_.internalField() - D_.prevIter().internalField())
                   /(maxDD + SMALL)
                );

            if (lduMatrix::debug)
            {
                Info << "Relative residual = " << residual << endl;
            }
        }
    }
    while(++iCorr < nCorr);

    U_ = fvc::ddt(D_);

    // Calculate second Piola-Kirchhoff stress
    {
        volSymmTensorField E = symm(gradD_);

        if(nonLinear)
        {
            E += 0.5*symm(gradD_ & gradD_.T());
        }

        sigma_ = 2*mu_*E + I*(lambda_*tr(E));
    }

    Info << solverPerf.solverName() << ": Solving for " << D_.name()
        << ", Initial residula = " << initialResidual
        << ", Final residual = " << solverPerf.initialResidual()
        << ", No outer iterations = " << iCorr
        << ", Relative momentum residual = " << residual << endl;

    lduMatrix::debug = 1;

    return residual;
}

scalar unsTotalLagrangianSolid::residual() const
{
    scalar res = GREAT;

    // Calculate momentum residual
    {
        scalar maxDD =
            gMax
            (
                mag
                (
                    D_.internalField()
                  - D_.oldTime().internalField()
                )
            );

        res =
            gMax
            (
                mag(D_.internalField() - D_.prevIter().internalField())
               /(maxDD + SMALL)
            );
    }

    return res;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidSolvers
} // End namespace Foam

// ************************************************************************* //
