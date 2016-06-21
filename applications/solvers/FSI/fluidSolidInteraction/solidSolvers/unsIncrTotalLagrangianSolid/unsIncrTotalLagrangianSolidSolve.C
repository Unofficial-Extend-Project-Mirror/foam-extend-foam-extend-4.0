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

#include "unsIncrTotalLagrangianSolid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGradf.H"
#include "tractionDisplacementIncrementFvPatchVectorField.H"
#include "skewCorrectionVectors.H"
#include "multiMaterial.H"
#include "twoDPointCorrector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidSolvers
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volVectorField> unsIncrTotalLagrangianSolid::residual
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

    bool enforceLinear = false;
    solidProperties().set("enforceLinear", enforceLinear);

    fvVectorMatrix DDEqn
    (
        rho_*fvm::d2dt2(DD_)
      - fvm::laplacian(2*muf_ + lambdaf_, DD_, "laplacian(DDD,DD)")
      + fvc::laplacian(muf_ + lambdaf_, DD_, "laplacian(DDD,DD)")
     == fvc::div
        (
            mesh().Sf()
          & (
                muf_*gradDDf_.T() + lambdaf_*(I*tr(gradDDf_))
            )
        )
    );

    // Add damping if not zero
    if (K.value() > SMALL)
    {
        DDEqn += K*rho_*fvm::ddt(DD_);
    }

    if (nonLinear && !enforceLinear)
    {
        surfaceSymmTensorField DEf = symm(gradDDf_);
        DEf += 0.5*symm(gradDDf_ & gradDDf_.T());
        DEf += 0.5*symm(gradDDf_ & gradDf_.T());
        DEf += 0.5*symm(gradDf_ & gradDDf_.T());

        DSigmaf_ = 2*muf_*DEf + I*(lambdaf_*tr(DEf));

        if (rheology_.plasticityActive())
        {
            DSigmaf_ -= 2*muf_*fvc::interpolate(rheology_.DEpsilonP());
        }

        DDEqn -=
            fvc::div
            (
                muf_*(mesh().Sf() & (gradDDf_ & gradDf_.T()))
              + muf_*(mesh().Sf() & (gradDf_ & gradDDf_.T()))
              + muf_*(mesh().Sf() & (gradDDf_ & gradDDf_.T()))
              + 0.5*lambdaf_*tr(gradDDf_ & gradDf_.T())*mesh().Sf()
              + 0.5*lambdaf_*tr(gradDf_ & gradDDf_.T())*mesh().Sf()
              + 0.5*lambdaf_*tr(gradDDf_ & gradDDf_.T())*mesh().Sf()
            )
          + fvc::div(mesh().Sf() & (DSigmaf_ & gradDf_))
          + fvc::div(mesh().Sf() & ((sigmaf_ + DSigmaf_) & gradDDf_));
    }

    if (rheology_.plasticityActive())
    {
        DDEqn +=
            fvc::div
            (
                2*muf_
               *(
                    mesh().Sf()
                  & fvc::interpolate(rheology_.DEpsilonP())
                )
            );
    }

    if (interface().valid())
    {
        interface()->correct(DDEqn);
    }

    DDEqn += source;

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
            dimensionedVector("0", dimForce, vector::zero),
            zeroGradientFvPatchVectorField::typeName
        )
    );

    tResidual().internalField() = DDEqn.residual();

    tResidual().correctBoundaryConditions();

    return tResidual;
}


void unsIncrTotalLagrangianSolid::initialise
(
    const volVectorField& sol,
    bool consistentBoundaryField
)
{
    Info << "Initialise solid solver solution" << endl;

    DD_ = sol;

    if (interface().valid())
    {
        interface()->updateDisplacement(pointDD_);
        interface()->updateDisplacementGradient(gradDD_, gradDDf_);

//         interface()->volToPointInterpolate
//         (
//             DD_,
//             interface()->displacementIncrement(),
//             pointDD_
//         );
    }
    else
    {
        volToPoint_.interpolate(DD_, pointDD_);
        gradDD_ = fvc::grad(DD_, pointDD_);
        gradDDf_ = fvc::fGrad(DD_, pointDD_);
    }

    DD_.correctBoundaryConditions();

    if (interface().valid())
    {
        interface()->updateDisplacement(pointDD_);
        interface()->updateDisplacementGradient(gradDD_, gradDDf_);
//         interface()->volToPointInterpolate
//         (
//             DD_,
//             interface()->displacementIncrement(),
//             pointDD_
//         );
    }
    else
    {
        volToPoint_.interpolate(DD_, pointDD_);
        gradDD_ = fvc::grad(DD_, pointDD_);
        gradDDf_ = fvc::fGrad(DD_, pointDD_);
    }

    DD_.correctBoundaryConditions();
}


scalar unsIncrTotalLagrangianSolid::smooth
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

    volScalarField rho = rheology_.rho();
    volScalarField mu = rheology_.mu();
    volScalarField lambda = rheology_.lambda();

    surfaceScalarField muf("muf", fvc::interpolate(rheology_.mu()));
    surfaceScalarField lambdaf
    (
        "lambdaf",
        fvc::interpolate(rheology_.lambda())
    );

    if (interface().valid())
    {
        muf_ = interface()->interpolate(mu);
        lambdaf_ = interface()->interpolate(lambda);

//         interface()->modifyProperty(muf);
//         interface()->modifyProperty(lambdaf);
    }

    int iCorr = 0;
    scalar initialResidual = 0;
    lduSolverPerformance solverPerf;
    scalar residual = GREAT;

    lduMatrix::debug = debug;

    bool enforceLinear = false;
    solidProperties().set("enforceLinear", enforceLinear);

    fvc::ddt(D_);

    do
    {
        if (lduMatrix::debug)
        {
            Info<< "Time: " << runTime().timeName()
                << ", outer iteration: " << iCorr << endl;
        }

        DD_.storePrevIter();

        fvVectorMatrix DDEqn
        (
            rho*fvm::d2dt2(DD_)
          - fvm::laplacian(2*muf + lambdaf, DD_, "laplacian(DDD,DD)")
          + fvc::laplacian(muf + lambdaf, DD_, "laplacian(DDD,DD)")
         == fvc::div
            (
                mesh().Sf()
              & (
                    muf*gradDDf_.T() + lambdaf*(I*tr(gradDDf_))
                )
            )
        );

        // Add damping if not zero
        if (K.value() > SMALL)
        {
            DDEqn += K*rho*fvm::ddt(DD_);
        }

        if (nonLinear && !enforceLinear)
        {
            surfaceSymmTensorField DEf = symm(gradDDf_);
            DEf += 0.5*symm(gradDDf_ & gradDDf_.T());
            DEf += 0.5*symm(gradDDf_ & gradDf_.T());
            DEf += 0.5*symm(gradDf_ & gradDDf_.T());

            DSigmaf_ = 2*muf*DEf + I*(lambdaf*tr(DEf));

            if (rheology_.plasticityActive())
            {
                DSigmaf_ -= 2*muf*fvc::interpolate(rheology_.DEpsilonP());
            }

            DDEqn -=
                fvc::div
                (
                    muf*(mesh().Sf() & (gradDDf_ & gradDf_.T()))
                  + muf*(mesh().Sf() & (gradDf_ & gradDDf_.T()))
                  + muf*(mesh().Sf() & (gradDDf_ & gradDDf_.T()))
                  + 0.5*lambdaf*tr(gradDDf_ & gradDf_.T())*mesh().Sf()
                  + 0.5*lambdaf*tr(gradDf_ & gradDDf_.T())*mesh().Sf()
                  + 0.5*lambdaf*tr(gradDDf_ & gradDDf_.T())*mesh().Sf()
                )
              + fvc::div(mesh().Sf() & (DSigmaf_ & gradDf_))
              + fvc::div(mesh().Sf() & ((sigmaf_ + DSigmaf_) & gradDDf_));
        }

        if (rheology_.plasticityActive())
        {
            DDEqn +=
                fvc::div
                (
                    2*muf
                   *(
                        mesh().Sf()
                      & fvc::interpolate(rheology_.DEpsilonP())
                    )
                );
        }

        if (interface().valid())
        {
            interface()->correct(DDEqn);
        }

        solverPerf = DDEqn.solve();

        if(iCorr == 0)
        {
            initialResidual = solverPerf.initialResidual();
        }

        DD_.relax();

        if (interface().valid())
        {
            interface()->updateDisplacement(pointDD_);
            interface()->updateDisplacementGradient(gradDD_, gradDDf_);

//             interface()->volToPointInterpolate
//             (
//                 DD_,
//                 interface()->displacementIncrement(),
//                 pointDD_
//             );
        }
        else
        {
            volToPoint_.interpolate(DD_, pointDD_);
            gradDD_ = fvc::grad(DD_, pointDD_);
            gradDDf_ = fvc::fGrad(DD_, pointDD_);
        }

        if (nonLinear && !enforceLinear)
        {
            surfaceScalarField Det = det(I+gradDf_+gradDDf_);

            scalar minDetFf = min(Det).value();
            reduce(minDetFf, minOp<scalar>());

            scalar maxDetFf = max(Det).value();
            reduce(maxDetFf, maxOp<scalar>());

            if ( (minDetFf<0.01) || (maxDetFf>100) )
            {
                enforceLinear = true;
                solidProperties().set("enforceLinear", enforceLinear);
            }
        }

        // Calculate momentu residual
        {
            scalar maxDD =
                gMax
                (
                    mag(DD_.internalField())
                );

            residual =
                gMax
                (
                    mag(DD_.internalField() - DD_.prevIter().internalField())
                   /(maxDD + SMALL)
                );

            if (lduMatrix::debug)
            {
                Info << "Relative residual = " << residual << endl;
            }
        }

        // Calculate strain increment
        {
            DEpsilon_ = symm(gradDD_);

            if(nonLinear && !enforceLinear)
            {
                DEpsilon_ += 0.5*symm(gradDD_ & gradDD_.T());
                DEpsilon_ += 0.5*symm(gradDD_ & gradD_.T());
                DEpsilon_ += 0.5*symm(gradD_ & gradDD_.T());
            }
        }

        // Correct plasticity term
        rheology_.correct();
    }
    while (++iCorr < nCorr);

//     DU_ = fvc::ddt(DD_);

    // Calculate second Piola-Kirchhoff stress increment
    {
        DSigma_ = 2*mu*DEpsilon_ + I*(lambda*tr(DEpsilon_));

        if (rheology_.plasticityActive())
        {
            DSigma_ -= 2*mu*rheology_.DEpsilonP();
        }
    }

    Info << solverPerf.solverName() << ": Solving for " << DD_.name()
        << ", Initial residula = " << initialResidual
        << ", Final residual = " << solverPerf.initialResidual()
        << ", No outer iterations = " << iCorr
        << ", Relative momentum residual = " << residual
        << ", enforceLinear = " << enforceLinear << endl;

    return residual;
}

scalar unsIncrTotalLagrangianSolid::smooth
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

    volScalarField rho = rheology_.rho();
    volScalarField mu = rheology_.mu();
    volScalarField lambda = rheology_.lambda();

    surfaceScalarField muf("muf", fvc::interpolate(rheology_.mu()));
    surfaceScalarField lambdaf
    (
        "lambdaf",
        fvc::interpolate(rheology_.lambda())
    );

    if (interface().valid())
    {
        muf_ = interface()->interpolate(mu);
        lambdaf_ = interface()->interpolate(lambda);

//         interface()->modifyProperty(muf);
//         interface()->modifyProperty(lambdaf);
    }

    int iCorr = 0;
    scalar initialResidual = 0;
    lduSolverPerformance solverPerf;
    scalar residual = GREAT;

    lduMatrix::debug = debug;

    bool enforceLinear = false;
    solidProperties().set("enforceLinear", enforceLinear);

    fvc::ddt(D_);

    do
    {
        if (lduMatrix::debug)
        {
            Info<< "Time: " << runTime().timeName()
                << ", outer iteration: " << iCorr << endl;
        }

        DD_.storePrevIter();

        fvVectorMatrix DDEqn
        (
            rho*fvm::d2dt2(DD_)
          - fvm::laplacian(2*muf + lambdaf, DD_, "laplacian(DDD,DD)")
          + fvc::laplacian(muf + lambdaf, DD_, "laplacian(DDD,DD)")
         == fvc::div
            (
                mesh().Sf()
              & (
                    muf*gradDDf_.T() + lambdaf*(I*tr(gradDDf_))
                )
            )
        );

        // Add damping if not zero
        if (K.value() > SMALL)
        {
            DDEqn += K*rho*fvm::ddt(DD_);
        }

        if (nonLinear && !enforceLinear)
        {
            surfaceSymmTensorField DEf = symm(gradDDf_);
            DEf += 0.5*symm(gradDDf_ & gradDDf_.T());
            DEf += 0.5*symm(gradDDf_ & gradDf_.T());
            DEf += 0.5*symm(gradDf_ & gradDDf_.T());

            DSigmaf_ = 2*muf*DEf + I*(lambdaf*tr(DEf));

            if (rheology_.plasticityActive())
            {
                DSigmaf_ -= 2*muf*fvc::interpolate(rheology_.DEpsilonP());
            }

            DDEqn -=
                fvc::div
                (
                    muf*(mesh().Sf() & (gradDDf_ & gradDf_.T()))
                  + muf*(mesh().Sf() & (gradDf_ & gradDDf_.T()))
                  + muf*(mesh().Sf() & (gradDDf_ & gradDDf_.T()))
                  + 0.5*lambdaf*tr(gradDDf_ & gradDf_.T())*mesh().Sf()
                  + 0.5*lambdaf*tr(gradDf_ & gradDDf_.T())*mesh().Sf()
                  + 0.5*lambdaf*tr(gradDDf_ & gradDDf_.T())*mesh().Sf()
                )
              + fvc::div(mesh().Sf() & (DSigmaf_ & gradDf_))
              + fvc::div(mesh().Sf() & ((sigmaf_ + DSigmaf_) & gradDDf_));
        }

        if (rheology_.plasticityActive())
        {
            DDEqn +=
                fvc::div
                (
                    2*muf
                   *(
                        mesh().Sf()
                      & fvc::interpolate(rheology_.DEpsilonP())
                    )
                );
        }

        if (interface().valid())
        {
            interface()->correct(DDEqn);
        }

        solverPerf = DDEqn.solve();

        if(iCorr == 0)
        {
            initialResidual = solverPerf.initialResidual();
        }

        DD_.relax();

        if (interface().valid())
        {
            interface()->updateDisplacement(pointDD_);
            interface()->updateDisplacementGradient
            (
                gradDD_,
                gradDDf_
            );

//             interface()->volToPointInterpolate
//             (
//                 DD_,
//                 interface()->displacementIncrement(),
//                 pointDD_
//             );
        }
        else
        {
            volToPoint_.interpolate(DD_, pointDD_);
            gradDD_ = fvc::grad(DD_, pointDD_);
            gradDDf_ = fvc::fGrad(DD_, pointDD_);
        }

        if (nonLinear && !enforceLinear)
        {
            surfaceScalarField Det = det(I+gradDf_+gradDDf_);

            scalar minDetFf = min(Det).value();
            reduce(minDetFf, minOp<scalar>());

            scalar maxDetFf = max(Det).value();
            reduce(maxDetFf, maxOp<scalar>());

            if ( (minDetFf<0.01) || (maxDetFf>100) )
            {
                enforceLinear = true;
                solidProperties().set("enforceLinear", enforceLinear);
            }
        }

        // Calculate momentu residual
        {
            scalar maxDD =
                gMax
                (
                    mag(DD_.internalField())
                );

            residual =
                gMax
                (
                    mag(DD_.internalField() - DD_.prevIter().internalField())
                   /(maxDD + SMALL)
                );

            if (lduMatrix::debug)
            {
                Info << "Relative residual = " << residual << endl;
            }
        }

        // Calculate strain increment
        {
            DEpsilon_ = symm(gradDD_);

            if(nonLinear && !enforceLinear)
            {
                DEpsilon_ += 0.5*symm(gradDD_ & gradDD_.T());
                DEpsilon_ += 0.5*symm(gradDD_ & gradD_.T());
                DEpsilon_ += 0.5*symm(gradD_ & gradDD_.T());
            }
        }

        // Correct plasticity term
        rheology_.correct();
    }
    while (++iCorr < nCorr);

//     DU_ = fvc::ddt(DD_);

    // Calculate second Piola-Kirchhoff stress increment
    {
        DSigma_ = 2*mu*DEpsilon_ + I*(lambda*tr(DEpsilon_));

        if (rheology_.plasticityActive())
        {
            DSigma_ -= 2*mu*rheology_.DEpsilonP();
        }
    }

    Info << solverPerf.solverName() << ": Solving for " << DD_.name()
        << ", Initial residula = " << initialResidual
        << ", Final residual = " << solverPerf.initialResidual()
        << ", No outer iterations = " << iCorr
        << ", Relative momentum residual = " << residual
        << ", enforceLinear = " << enforceLinear << endl;

    return residual;
}


scalar unsIncrTotalLagrangianSolid::residual() const
{
    scalar res = GREAT;

    // Calculate momentum residual
    {
        scalar maxDD =
            gMax
            (
                mag(DD_.internalField())
            );

        res =
            gMax
            (
                mag(DD_.internalField() - DD_.prevIter().internalField())
               /(maxDD + SMALL)
            );
    }

    return res;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidSolvers
} // End namespace Foam

// ************************************************************************* //
