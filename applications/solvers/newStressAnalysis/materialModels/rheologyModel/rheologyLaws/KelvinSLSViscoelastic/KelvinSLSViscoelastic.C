/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description

\*---------------------------------------------------------------------------*/

#include "KelvinSLSViscoelastic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(KelvinSLSViscoelastic, 0);
    addToRunTimeSelectionTable(rheologyLaw, KelvinSLSViscoelastic, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::KelvinSLSViscoelastic::KelvinSLSViscoelastic
(
    const word& name,
    const volSymmTensorField& sigma,
    const dictionary& dict
)
:
    rheologyLaw(name, sigma, dict),
    rho_(dict.lookup("rho")),
    k1_(dict.lookup("k1")),
    k2_(dict.lookup("k2")),
    eta2_(dict.lookup("eta2")),
    nu_(dict.lookup("nu"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::KelvinSLSViscoelastic::~KelvinSLSViscoelastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::KelvinSLSViscoelastic::rho(scalar t) const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::KelvinSLSViscoelastic::E(scalar t) const
{
    scalar E = 0.0;

    if(t>=0)
    {
        scalar p1 = eta2_.value()/(k1_.value() + k2_.value());
 
        scalar q0 = k1_.value()*k2_.value()/(k1_.value() + k2_.value());

        scalar q1 = k1_.value()*eta2_.value()/(k1_.value() + k2_.value());

        E = q0 + (q1/p1 - q0)*exp(-t/p1);
    }
    

    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("E", k1_.dimensions(), E),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::KelvinSLSViscoelastic::nu(scalar t) const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "nu",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            nu_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::KelvinSLSViscoelastic::J(scalar t) const
{
    scalar J = 0.0;

    if(t >= 0)
    {
        scalar p1 = eta2_.value()/(k1_.value() + k2_.value());
 
        scalar q0 = k1_.value()*k2_.value()/(k1_.value() + k2_.value());

        scalar q1 = k1_.value()*eta2_.value()/(k1_.value() + k2_.value());

        J = 1.0/q0 + (p1/q1 - 1.0/q0)*exp(-q0*t/q1);
    }

    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "J",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("J", dimless/k1_.dimensions(), J),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();

    return tresult;
}


// ************************************************************************* //
