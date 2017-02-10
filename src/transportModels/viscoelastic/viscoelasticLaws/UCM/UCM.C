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

#include "UCM.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(UCM, 0);
    addToRunTimeSelectionTable(viscoelasticLaw, UCM, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UCM::UCM
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
    etaP_(dict.lookup("etaP")),
    lambda_(dict.lookup("lambda")),
    etaStab_(dimensionedScalar::lookupOrDefault("etaStab", dict, 0.0, dimMass/(dimLength*dimTime)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::UCM::divTau(volVectorField& U) const
{
    if(etaStab_.value() < SMALL)
    {
        dimensionedScalar etaPEff = etaP_;

    	return
    	(
    	    fvc::div(tau_/rho_, "div(tau)")
    	  - fvc::laplacian(etaPEff/rho_, U, "laplacian(etaPEff,U)")
    	  + fvm::laplacian( (etaPEff)/rho_, U, "laplacian(etaPEff+etaS,U)")
    	);
    }
    else
    {
    	return
    	(
    	    fvc::div(tau_/rho_, "div(tau)")
    	  - fvc::div((etaStab_/rho_)*fvc::grad(U), "div(tau)")
    	  + fvm::laplacian( (etaStab_)/rho_, U, "laplacian(etaPEff+etaS,U)")
    	);
    }
}


void Foam::UCM::correct()
{
    // Velocity gradient tensor
    volTensorField L = fvc::grad(U());

    // Convected derivate term
    volTensorField C = tau_ & L;

    // Twice the rate of deformation tensor
    volSymmTensorField twoD = twoSymm(L);

     // Stress transport equation
    fvSymmTensorMatrix tauEqn
    (
        fvm::ddt(tau_)
      + fvm::div(phi(), tau_)
     ==
        etaP_/lambda_*twoD
      + twoSymm(C)
      - fvm::Sp( 1/lambda_, tau_)
    );

    tauEqn.relax();
    tauEqn.solve();
}


// ************************************************************************* //
