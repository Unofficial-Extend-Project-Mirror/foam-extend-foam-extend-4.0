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

#include "Feta_PTT.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Feta_PTT, 0);
    addToRunTimeSelectionTable(viscoelasticLaw, Feta_PTT, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Feta_PTT::Feta_PTT
(
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    viscoelasticLaw(name, U, phi),
    tau_
    (
        IOobject
        (
            "tau" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    rho_(dict.lookup("rho")),
    etaS_(dict.lookup("etaS")),
    etaP_(dict.lookup("etaP")),
    epsilon_(dict.lookup("epsilon")),
    lambda_(dict.lookup("lambda")),
    zeta_(dict.lookup("zeta")),
    A_(dict.lookup("A")),
    a_(dict.lookup("a")),
    b_(dict.lookup("b")),
    etaPEff_
    (
        IOobject
        (
            "etaPEff",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        etaP_/
        (
            Foam::pow(scalar(1) + A_*Foam::pow(0.5*( Foam::sqr(tr(tau_))
          - tr(tau_ & tau_))*Foam::sqr(lambda_)/Foam::sqr(etaP_), a_), b_)
        )
    ),
    lambdaEff_
    (
        IOobject
        (
            "lambdaEff",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        lambda_ / (scalar(1)  + epsilon_*lambda_*tr(tau_) / etaP_)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::Feta_PTT::divTau(volVectorField& U) const
{
    // Need to be equal to old time step (a constant)
    dimensionedScalar etaPEff = etaP_;

    return
    (
        fvc::div(tau_/rho_, "div(tau)")
      - fvc::laplacian(etaPEff/rho_, U, "laplacian(etaPEff,U)")
      + fvm::laplacian( (etaPEff + etaS_)/rho_, U, "laplacian(etaPEff+etaS,U)")
    );
}


void Foam::Feta_PTT::correct()
{
    // Velocity gradient tensor
    volTensorField L = fvc::grad(U());

    // Convected derivate term
    volTensorField C = tau_ & L;

    // Twice the rate of deformation tensor
    volSymmTensorField twoD = twoSymm(L);

    // etaP effective
    etaPEff_ = etaP_/
        (
            Foam::pow(scalar(1) + A_*Foam::pow(0.5*( Foam::sqr(tr(tau_))
          - tr(tau_ & tau_)) * Foam::sqr(lambda_)/Foam::sqr(etaP_), a_), b_)
        );

    // lambda effective
    lambdaEff_ = lambda_/(scalar(1)  + epsilon_*lambda_*tr(tau_)/etaP_);

    // Stress transport equation
    fvSymmTensorMatrix tauEqn
    (
        fvm::ddt(tau_)
      + fvm::div(phi(), tau_)
     ==
        etaPEff_/lambdaEff_*twoD
      + twoSymm(C)
      - zeta_/2*((tau_ & twoD) + (twoD & tau_))
      - fvm::Sp(epsilon_/etaPEff_*tr(tau_) + 1/lambdaEff_, tau_)
    );

    tauEqn.relax();
    tauEqn.solve();
}


// ************************************************************************* //
