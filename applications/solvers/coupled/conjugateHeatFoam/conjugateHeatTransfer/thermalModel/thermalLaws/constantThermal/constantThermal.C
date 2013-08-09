/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 Hrvoje Jasak
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

\*---------------------------------------------------------------------------*/

#include "constantThermal.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantThermal, 0);
    addToRunTimeSelectionTable(thermalLaw, constantThermal, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::constantThermal::constantThermal
(
    const word& name,
    const volScalarField& T,
    const dictionary& dict
)
:
    thermalLaw(name, T, dict),
    rho_(dict.lookup("rho")),
    C_(dict.lookup("C")),
    k_(dict.lookup("k"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantThermal::~constantThermal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::constantThermal::rho() const
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

Foam::tmp<Foam::volScalarField> Foam::constantThermal::C() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "C",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            C_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::constantThermal::k() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "kTmp",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            k_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tresult().correctBoundaryConditions();

    return tresult;
}


// ************************************************************************* //
