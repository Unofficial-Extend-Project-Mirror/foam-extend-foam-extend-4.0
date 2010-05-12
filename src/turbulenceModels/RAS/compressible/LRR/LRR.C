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

#include "LRR.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LRR, 0);
addToRunTimeSelectionTable(RASModel, LRR, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from components
LRR::LRR
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    basicThermo& thermophysicalModel
)
:
    RASModel(typeName, rho, U, phi, thermophysicalModel),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    Clrr1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clrr1",
            coeffDict_,
            1.8
        )
    ),
    Clrr2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clrr2",
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
    couplingFactor_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "couplingFactor",
            coeffDict_,
            0.0
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
    alphah_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphah",
            coeffDict_,
            1.0
        )
    ),

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

    mut_
    (
        IOobject
        (
            "mut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Cmu_*rho_*sqr(k_)/(epsilon_ + epsilonSmall_)
    )
{
#   include "wallViscosityI.H"

    if (couplingFactor_.value() < 0.0 || couplingFactor_.value() > 1.0)
    {
        FatalErrorIn
        (
            "LRR::LRR"
            "( const volScalarField&, const volVectorField&"
            ", const surfaceScalarField&, incompressibleTransportModel&)"
        )   << "couplingFactor = " << couplingFactor_
            << " is not in range 0 - 1" << nl
            << exit(FatalError);
    }

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> LRR::devRhoReff() const
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
            rho_*R_ - mu()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> LRR::divDevRhoReff(volVectorField& U) const
{
    if (couplingFactor_.value() > 0.0)
    {
        return
        (
            fvc::div(rho_*R_ + couplingFactor_*mut_*fvc::grad(U))
          + fvc::laplacian((1.0 - couplingFactor_)*mut_, U)
          - fvm::laplacian(muEff(), U)
          - fvc::div(mu()*dev2(fvc::grad(U)().T()))
        );
    }
    else
    {
        return
        (
            fvc::div(rho_*R_)
          + fvc::laplacian(mut_, U)
          - fvm::laplacian(muEff(), U)
          - fvc::div(mu()*dev2(fvc::grad(U)().T()))
        );
    }
}


bool LRR::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict_);
        Clrr1_.readIfPresent(coeffDict_);
        Clrr2_.readIfPresent(coeffDict_);
        C1_.readIfPresent(coeffDict_);
        C2_.readIfPresent(coeffDict_);
        Cs_.readIfPresent(coeffDict_);
        Ceps_.readIfPresent(coeffDict_);
        alphaR_.readIfPresent(coeffDict_);
        alphaEps_.readIfPresent(coeffDict_);
        alphah_.readIfPresent(coeffDict_);
        couplingFactor_.readIfPresent(coeffDict_);

        if (couplingFactor_.value() < 0.0 || couplingFactor_.value() > 1.0)
        {
            FatalErrorIn("LRR::read()")
                << "couplingFactor = " << couplingFactor_
                << " is not in range 0 - 1" << nl
                << exit(FatalError);
        }

        return true;
    }
    else
    {
        return false;
    }
}


void LRR::correct()
{
    if (!turbulence_)
    {
        // Re-calculate viscosity
        mut_ = rho_*Cmu_*sqr(k_)/(epsilon_ + epsilonSmall_);
        return;
    }

    RASModel::correct();

    volSymmTensorField P = -twoSymm(R_ & fvc::grad(U_));
    volScalarField G = 0.5*mag(tr(P));

#   include "wallFunctionsI.H"

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(rho_, epsilon_)
      + fvm::div(phi_, epsilon_)
    //- fvm::laplacian(Ceps*rho_*(k_/epsilon_)*R_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        C1_*rho_*G*epsilon_/k_
      - fvm::Sp(C2_*rho_*epsilon_/k_, epsilon_)
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
                P[faceCelli]
                    *= min(G[faceCelli]/(0.5*mag(tr(P[faceCelli])) + SMALL), 100.0);
            }
        }
    }


    tmp<fvSymmTensorMatrix> REqn
    (
        fvm::ddt(rho_, R_)
      + fvm::div(phi_, R_)
    //- fvm::laplacian(Cs*rho_*(k_/epsilon_)*R_, R_)
      - fvm::laplacian(DREff(), R_)
      + fvm::Sp(Clrr1_*rho_*epsilon_/k_, R_)
     ==
        rho_*P
      - (2.0/3.0*(1 - Clrr1_)*I)*rho_*epsilon_
      - Clrr2_*rho_*dev(P)
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

    k_ = 0.5*tr(R_);
    bound(k_, k0_);


    // Re-calculate viscosity
    mut_ = rho_*Cmu_*sqr(k_)/epsilon_;

#   include "wallViscosityI.H"


    // Correct wall shear stresses

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (typeid(curPatch) == typeid(wallFvPatch))
        {
            symmTensorField& Rw = R_.boundaryField()[patchi];

            const scalarField& rhow = rho_.boundaryField()[patchi];
            const scalarField& mutw = mut_.boundaryField()[patchi];

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
                tensor tauw = -(mutw[facei]/rhow[facei])*2*dev(symm(gradUw));

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
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
