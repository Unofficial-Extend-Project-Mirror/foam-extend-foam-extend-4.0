/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "LaunderGibsonRSTM.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LaunderGibsonRSTM, 0);
addToRunTimeSelectionTable(RASModel, LaunderGibsonRSTM, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
LaunderGibsonRSTM::LaunderGibsonRSTM
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    Clg1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clg1",
            coeffDict_,
            1.8
        )
    ),
    Clg2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clg2",
            coeffDict_,
            0.6
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.92
        )
    ),
    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            coeffDict_,
            0.25
        )
    ),
    Ceps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps",
            coeffDict_,
            0.15
        )
    ),
    alphaR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaR",
            coeffDict_,
            1.22
        )
    ),
    alphaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaEps",
            coeffDict_,
            0.76923
        )
    ),
    C1Ref_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1Ref",
            coeffDict_,
            0.5
        )
    ),
    C2Ref_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2Ref",
            coeffDict_,
            0.3
        )
    ),
    couplingFactor_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "couplingFactor",
            coeffDict_,
            0.0
        )
    ),

    yr_(mesh_),

    R_
    (
        IOobject
        (
            "R",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    nut_(Cmu_*sqr(k_)/(epsilon_ + epsilonSmall_))
{
#   include "wallViscosityI.H"

    if (couplingFactor_.value() < 0.0 || couplingFactor_.value() > 1.0)
    {
        FatalErrorIn
        (
            "LaunderGibsonRSTM::LaunderGibsonRSTM"
            "(const volVectorField& U, const surfaceScalarField& phi,"
            "transportModel& lamTransportModel)"
        )   << "couplingFactor = " << couplingFactor_
            << " is not in range 0 - 1" << nl
            << exit(FatalError);
    }

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> LaunderGibsonRSTM::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            R_ - nu()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> LaunderGibsonRSTM::divDevReff(volVectorField& U) const
{
    if (couplingFactor_.value() > 0.0)
    {
        return
        (
            fvc::div(R_ + couplingFactor_*nut_*fvc::grad(U), "div(R)")
          + fvc::laplacian((1.0-couplingFactor_)*nut_, U, "laplacian(nuEff,U)")
          - fvm::laplacian(nuEff(), U)
        );
    }
    else
    {
        return
        (
            fvc::div(R_)
          + fvc::laplacian(nut_, U, "laplacian(nuEff,U)")
          - fvm::laplacian(nuEff(), U)
        );
    }
}


bool LaunderGibsonRSTM::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict_);
        Clg1_.readIfPresent(coeffDict_);
        Clg2_.readIfPresent(coeffDict_);
        C1_.readIfPresent(coeffDict_);
        C2_.readIfPresent(coeffDict_);
        Cs_.readIfPresent(coeffDict_);
        Ceps_.readIfPresent(coeffDict_);
        alphaR_.readIfPresent(coeffDict_);
        alphaEps_.readIfPresent(coeffDict_);
        C1Ref_.readIfPresent(coeffDict_);
        C2Ref_.readIfPresent(coeffDict_);

        couplingFactor_.readIfPresent(coeffDict_);

        if (couplingFactor_.value() < 0.0 || couplingFactor_.value() > 1.0)
        {
            FatalErrorIn("LaunderGibsonRSTM::read()")
                << "couplingFactor = " << couplingFactor_
                << " is not in range 0 - 1"
                << exit(FatalError);
        }

        return true;
    }
    else
    {
        return false;
    }
}


void LaunderGibsonRSTM::correct()
{
    transportModel_.correct();

    if (!turbulence_)
    {
        return;
    }

    RASModel::correct();

    if (mesh_.changing())
    {
        yr_.correct();
    }

    volSymmTensorField P = -twoSymm(R_ & fvc::grad(U_));
    volScalarField G = 0.5*mag(tr(P));

#   include "wallFunctionsI.H"

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      + fvm::SuSp(-fvc::div(phi_), epsilon_)
    //- fvm::laplacian(Ceps*(k_/epsilon_)*R_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        C1_*G*epsilon_/k_
      - fvm::Sp(C2_*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();

#   include "wallDissipationI.H"

    solve(epsEqn);
    bound(epsilon_, epsilon0_);


    // Reynolds stress equation

    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (typeid(curPatch) == typeid(wallFvPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];
                P[faceCelli] *=
                    min(G[faceCelli]/(0.5*mag(tr(P[faceCelli])) + SMALL), 1.0);
            }
        }
    }

    volSymmTensorField reflect = C1Ref_*epsilon_/k_*R_ - C2Ref_*Clg2_*dev(P);

    tmp<fvSymmTensorMatrix> REqn
    (
        fvm::ddt(R_)
      + fvm::div(phi_, R_)
      + fvm::SuSp(-fvc::div(phi_), R_)
    //- fvm::laplacian(Cs*(k_/epsilon_)*R_, R_)
      - fvm::laplacian(DREff(), R_)
      + fvm::Sp(Clg1_*epsilon_/k_, R_)
      ==
        P
      + (2.0/3.0*(Clg1_ - 1)*I)*epsilon_
      - Clg2_*dev(P)

        // wall reflection terms
      + symm
        (
            I*((yr_.n() & reflect) & yr_.n())
          - 1.5*(yr_.n()*(reflect & yr_.n())
          + (yr_.n() & reflect)*yr_.n())
        )*pow(Cmu_, 0.75)*pow(k_, 1.5)/(kappa_*yr_*epsilon_)
    );

    REqn().relax();
    solve(REqn);

    R_.max
    (
        dimensionedSymmTensor
        (
            "zero",
            R_.dimensions(),
            symmTensor
            (
                k0_.value(), -GREAT, -GREAT,
                             k0_.value(), -GREAT,
                                          k0_.value()
            )
        )
    );

    k_ == 0.5*tr(R_);
    bound(k_, k0_);


    // Re-calculate turbulent viscosity
    nut_ = Cmu_*sqr(k_)/epsilon_;


#   include "wallViscosityI.H"


    // Correct wall shear stresses

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (typeid(curPatch) == typeid(wallFvPatch))
        {
            symmTensorField& Rw = R_.boundaryField()[patchi];

            const scalarField& nutw = nut_.boundaryField()[patchi];

            vectorField snGradU = U_.boundaryField()[patchi].snGrad();

            const vectorField& faceAreas
                = mesh_.Sf().boundaryField()[patchi];

            const scalarField& magFaceAreas
                = mesh_.magSf().boundaryField()[patchi];

            forAll(curPatch, facei)
            {
                // Calculate near-wall velocity gradient
                tensor gradUw
                    = (faceAreas[facei]/magFaceAreas[facei])*snGradU[facei];

                // Calculate near-wall shear-stress tensor
                tensor tauw = -nutw[facei]*2*symm(gradUw);

                // Reset the shear components of the stress tensor
                Rw[facei].xy() = tauw.xy();
                Rw[facei].xz() = tauw.xz();
                Rw[facei].yz() = tauw.yz();
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
