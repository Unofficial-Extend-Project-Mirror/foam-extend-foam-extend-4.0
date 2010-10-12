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

#include "S_MDCPP.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(S_MDCPP, 0);
    addToRunTimeSelectionTable(viscoelasticLaw, S_MDCPP, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::S_MDCPP::S_MDCPP
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

Foam::tmp<Foam::fvVectorMatrix> Foam::S_MDCPP::divTau(volVectorField& U) const
{
    dimensionedScalar etaPEff = etaP_;

    return
    (
        fvc::div(tau_/rho_, "div(tau)")
      - fvc::laplacian(etaPEff/rho_, U, "laplacian(etaPEff,U)")
      + fvm::laplacian( (etaPEff + etaS_)/rho_, U, "laplacian(etaPEff+etaS,U)")
    );

}


void Foam::S_MDCPP::correct()
{
    // Velocity gradient tensor
    volTensorField L = fvc::grad(U());

    // Convected derivate term
    volTensorField C = tau_ & L;

    // Twice the rate of deformation tensor
    volSymmTensorField twoD = twoSymm(L);

    // Lambda (Backbone stretch)
    volScalarField Lambda =
        Foam::sqrt(1 + tr(tau_)*lambdaOb_*(1 - zeta_)/3/etaP_);

    // Auxiliary field
    volScalarField aux = Foam::exp( 2/q_*(Lambda - 1));

    // Extra function
    volScalarField fTau =
        aux*(2*lambdaOb_/lambdaOs_*(1 - 1/Lambda) + 1/Foam::sqr(Lambda));

    // Stress transport equation
    fvSymmTensorMatrix tauEqn
    (
        fvm::ddt(tau_)
      + fvm::div(phi(), tau_)
     ==
        etaP_/lambdaOb_*twoD
      + twoSymm(C)
      - zeta_/2*((tau_ & twoD) + (twoD & tau_))
      - fvm::Sp(1/lambdaOb_*fTau, tau_)
      - (
            1/lambdaOb_*(etaP_/lambdaOb_/(1 - zeta_)*(fTau - aux)*I_)
        )
    );

    tauEqn.relax();
    tauEqn.solve();
}


// ************************************************************************* //
