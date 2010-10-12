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

#include "XPP_SE.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(XPP_SE, 0);
    addToRunTimeSelectionTable(viscoelasticLaw, XPP_SE, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XPP_SE::XPP_SE
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
    alpha_(dict.lookup("alpha")),
    lambdaOb_(dict.lookup("lambdaOb")),
    lambdaOs_(dict.lookup("lambdaOs")),
    q_(dict.lookup("q"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::XPP_SE::divTau(volVectorField& U) const
{

     dimensionedScalar etaPEff = etaP_;

    return
    (
        fvc::div(tau_/rho_, "div(tau)")
      - fvc::laplacian(etaPEff/rho_, U, "laplacian(etaPEff,U)")
      + fvm::laplacian( (etaPEff + etaS_)/rho_, U, "laplacian(etaPEff+etaS,U)")
    );

}


void Foam::XPP_SE::correct()
{
    // Velocity gradient tensor
    volTensorField L = fvc::grad(U());

    // Convected derivate term
    volTensorField C = tau_ & L;

    // Twice the rate of deformation tensor
    volSymmTensorField twoD = twoSymm(L);

    // Lambda (Backbone stretch)
    volScalarField Lambda = Foam::sqrt(1 + tr(tau_)*lambdaOb_/3/etaP_);

    // lambdaS (stretch relaxation time)
    volScalarField lambdaS = lambdaOs_*Foam::exp(-2/q_*(Lambda - 1));

    // Extra function
    volScalarField fTau = 2*lambdaOb_/lambdaS*(1 - 1/Lambda)
      + 1/Foam::sqr(Lambda)*
        (1 - alpha_*tr(tau_ & tau_)/3/Foam::sqr(etaP_/lambdaOb_));

     // Stress transport equation
    fvSymmTensorMatrix tauEqn
    (
        fvm::ddt(tau_)
      + fvm::div(phi(), tau_)
     ==
        etaP_/lambdaOb_*twoD
      + twoSymm(C)
      - fvm::Sp(1/lambdaOb_*fTau, tau_)
      - (
            1/lambdaOb_*
            (
                alpha_*lambdaOb_/etaP_*(tau_ & tau_)
              + etaP_/lambdaOb_*(fTau - 1)*I_
            )
        )
    );

    tauEqn.relax();
    tauEqn.solve();

}


// ************************************************************************* //
