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

#include "MaxwellElasticViscoelastic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MaxwellElasticViscoelastic, 0);
    addToRunTimeSelectionTable
    (
        rheologyLaw,
        MaxwellElasticViscoelastic,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::MaxwellElasticViscoelastic::MaxwellElasticViscoelastic
(
    const word& name,
    const volSymmTensorField& sigma,
    const dictionary& dict
)
:
    rheologyLaw(name, sigma, dict),
    rho_(dict.lookup("rho")),
    k_(dict.lookup("k")),
    eta_(dict.lookup("eta")),
    nu_(dict.lookup("nu"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::MaxwellElasticViscoelastic::~MaxwellElasticViscoelastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::MaxwellElasticViscoelastic::rho(scalar t) const
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


Foam::tmp<Foam::volScalarField>
Foam::MaxwellElasticViscoelastic::E(scalar t) const
{
    scalar tau = eta_.value()/k_.value();

    tmp<volScalarField> tE
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
            k_*exp(-t/tau),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (t < 0)
    {
        tE().internalField() = 0.0;
        tE().correctBoundaryConditions();
    }

    return tE;
}


Foam::tmp<Foam::volScalarField>
Foam::MaxwellElasticViscoelastic::nu(scalar t) const
{
    scalar tau = eta_.value()/k_.value();

    scalar E0 = k_.value();

    scalar Et = k_.value()*exp(-t/tau);

    scalar nu = 0.5 - Et/(6*E0);

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
            dimensionedScalar("nu", dimless, nu),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField>
Foam::MaxwellElasticViscoelastic::J(scalar t) const
{
    tmp<volScalarField> tJ
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
            dimensionedScalar
            (
                "J", 
                dimless/k_.dimensions(), 
                1.0/k_.value() + t/eta_.value()
            ),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (t < 0)
    {
        tJ().internalField() = 0.0;
        tJ().correctBoundaryConditions();
    }

    return tJ;
}


// ************************************************************************* //
