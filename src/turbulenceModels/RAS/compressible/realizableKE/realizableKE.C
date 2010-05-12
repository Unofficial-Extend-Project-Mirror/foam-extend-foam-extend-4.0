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

#include "realizableKE.H"
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

defineTypeNameAndDebug(realizableKE, 0);
addToRunTimeSelectionTable(RASModel, realizableKE, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> realizableKE::rCmu
(
    const volTensorField& gradU,
    const volScalarField& S2,
    const volScalarField& magS
)
{
    tmp<volSymmTensorField> tS = dev(symm(gradU));
    const volSymmTensorField& S = tS();

    volScalarField W =
        (2*sqrt(2.0))*((S&S)&&S)
       /(
            magS*S2
          + dimensionedScalar("small", dimensionSet(0, 0, -3, 0, 0), SMALL)
        );

    tS.clear();

    volScalarField phis =
        (1.0/3.0)*acos(min(max(sqrt(6.0)*W, -scalar(1)), scalar(1)));
    volScalarField As = sqrt(6.0)*cos(phis);
    volScalarField Us = sqrt(S2/2.0 + magSqr(skew(gradU)));

    return 1.0/(A0_ + As*Us*k_/(epsilon_ + epsilonSmall_));
}


tmp<volScalarField> realizableKE::rCmu
(
    const volTensorField& gradU
)
{
    volScalarField S2 = 2*magSqr(dev(symm(gradU)));
    volScalarField magS = sqrt(S2);
    return rCmu(gradU, S2, magS);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

realizableKE::realizableKE
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
    A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            coeffDict_,
            4.0
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.9
        )
    ),
    alphak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphak",
            coeffDict_,
            1.0
        )
    ),
    alphaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaEps",
            coeffDict_,
            0.833333
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
        rCmu(fvc::grad(U_))*rho_*sqr(k_)/(epsilon_ + epsilonSmall_)
    )
{
    bound(k_, k0_);
    bound(epsilon_, epsilon0_);
#   include "wallViscosityI.H"

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> realizableKE::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - (mut_/rho_)*dev(twoSymm(fvc::grad(U_))),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> realizableKE::devRhoReff() const
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
           -muEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> realizableKE::divDevRhoReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(muEff(), U) - fvc::div(muEff()*dev2(fvc::grad(U)().T()))
    );
}


bool realizableKE::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict_);
        A0_.readIfPresent(coeffDict_);
        C2_.readIfPresent(coeffDict_);
        alphak_.readIfPresent(coeffDict_);
        alphaEps_.readIfPresent(coeffDict_);
        alphah_.readIfPresent(coeffDict_);

        return true;
    }
    else
    {
        return false;
    }
}


void realizableKE::correct()
{
    if (!turbulence_)
    {
        // Re-calculate viscosity
        mut_ = rCmu(fvc::grad(U_))*rho_*sqr(k_)/epsilon_;
        return;
    }

    RASModel::correct();

    volScalarField divU = fvc::div(phi_/fvc::interpolate(rho_));

    if (mesh_.moving())
    {
        divU += fvc::div(mesh_.phi());
    }

    volTensorField gradU = fvc::grad(U_);
    volScalarField S2 = 2*magSqr(dev(symm(gradU)));
    volScalarField magS = sqrt(S2);

    volScalarField eta = magS*k_/epsilon_;
    volScalarField C1 = max(eta/(scalar(5) + eta), scalar(0.43));

    volScalarField G = mut_*(gradU && dev(twoSymm(gradU)));

#   include "wallFunctionsI.H"

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(rho_, epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        C1*rho_*magS*epsilon_
      - fvm::Sp
        (
            C2_*rho_*epsilon_/(k_ + sqrt((mu()/rho_)*epsilon_)),
            epsilon_
        )
    );

    epsEqn().relax();

#   include "wallDissipationI.H"

    solve(epsEqn);
    bound(epsilon_, epsilon0_);


    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(rho_, k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G - fvm::SuSp(2.0/3.0*rho_*divU, k_)
      - fvm::Sp(rho_*(epsilon_)/k_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);

    // Re-calculate viscosity
    mut_ = rCmu(gradU, S2, magS)*rho_*sqr(k_)/epsilon_;

#   include "wallViscosityI.H"

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
