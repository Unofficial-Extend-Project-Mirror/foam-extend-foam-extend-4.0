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

#include "DCPP.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DCPP, 0);
    addToRunTimeSelectionTable(viscoelasticLaw, DCPP, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DCPP::DCPP
(
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    viscoelasticLaw(name, U, phi),
    S_
    (
        IOobject
        (
            "S" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    Lambda_
    (
        IOobject
        (
            "Lambda" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    tau_
    (
        IOobject
        (
            "tau" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedSymmTensor
        (
            "zero",
            dimensionSet(1, -1, -2, 0, 0, 0, 0),
            symmTensor::zero
        )
    ),
    I_
    (
        dimensionedSymmTensor
        (
            "I",
            dimensionSet(0, 0, 0, 0, 0, 0, 0),
            symmTensor
            (
                1, 0, 0,
                   1, 0,
                      1
            )
        )
    ),
    rho_(dict.lookup("rho")),
    etaS_(dict.lookup("etaS")),
    etaP_(dict.lookup("etaP")),
    zeta_(dict.lookup("zeta")),
    lambdaOb_(dict.lookup("lambdaOb")),
    lambdaOs_(dict.lookup("lambdaOs")),
    q_(dict.lookup("q"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::DCPP::divTau(volVectorField& U) const
{
    dimensionedScalar etaPEff = etaP_;

    return
    (
        fvc::div(tau_/rho_, "div(tau)")
      - fvc::laplacian(etaPEff/rho_, U, "laplacian(etaPEff,U)")
      + fvm::laplacian((etaPEff + etaS_)/rho_, U, "laplacian(etaPEff+etaS,U)")
    );
}


void Foam::DCPP::correct()
{
    // Velocity gradient tensor
    volTensorField L = fvc::grad(U());

    // Upper convected derivate term
    volTensorField Cupper = S_ & L;

    // Lower convected derivate term
    volTensorField Clower = L & S_;

    // Twice the rate of deformation tensor
    volSymmTensorField twoD = twoSymm(L);


    // Evolution of orientation
    fvSymmTensorMatrix SEqn
    (
        fvm::ddt(S_)
      + fvm::div(phi(), S_)
     ==
        (1 - zeta_/2)*twoSymm(Cupper)
      - (zeta_/2)*twoSymm(Clower)
      - (1 - zeta_)*fvm::Sp((twoD && S_), S_)
      - fvm::Sp(1/lambdaOb_/Foam::sqr(Lambda_), S_)
      + 1/lambdaOb_/Foam::sqr(Lambda_)/3*I_
    );

    SEqn.relax();
    SEqn.solve();

     // Evolution of the backbone stretch
    fvScalarMatrix lambdaEqn
    (
        fvm::ddt(Lambda_)
      + fvm::div(phi(), Lambda_)
     ==
        fvm::Sp((twoD && S_)/2 , Lambda_)
      - fvm::Sp(Foam::exp(2/q_*(Lambda_ - 1))/lambdaOs_ , Lambda_)
      + Foam::exp( 2/q_*(Lambda_ - 1))/lambdaOs_
    );

    lambdaEqn.relax();
    lambdaEqn.solve();

    // Viscoelastic stress
    tau_ = etaP_/lambdaOb_/(1 - zeta_) * (3*Foam::sqr(Lambda_)*S_ - I_);
}


// ************************************************************************* //
