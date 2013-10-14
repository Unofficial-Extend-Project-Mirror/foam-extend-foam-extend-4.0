/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008 Hrvoje Jasak
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

#include "multiMaterialThermal.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiMaterialThermal, 0);
    addToRunTimeSelectionTable(thermalLaw, multiMaterialThermal, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::multiMaterialThermal::indicator
(
    const label i
) const
{
    const scalarField& mat = materials_.internalField();

    tmp<scalarField> tresult(new scalarField(mat.size(), 0.0));
    scalarField& result = tresult();

    forAll (mat, matI)
    {
        if (mat[matI] > i - SMALL && mat[matI] < i + 1 - SMALL)
        {
            result[matI] = 1.0;
        }
    }

    return tresult;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::multiMaterialThermal::multiMaterialThermal
(
    const word& name,
    const volScalarField& T,
    const dictionary& dict
)
:
    thermalLaw(name, T, dict),
    PtrList<thermalLaw>(),
    materials_
    (
        IOobject
        (
            "materials",
            mesh().time().timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    )
{
    PtrList<thermalLaw>& laws = *this;

    PtrList<entry> lawEntries(dict.lookup("laws"));
    laws.setSize(lawEntries.size());

    forAll (laws, lawI)
    {
        laws.set
        (
            lawI,
            thermalLaw::New
            (
                lawEntries[lawI].keyword(),
                T,
                lawEntries[lawI].dict()
            )
        );
    }

    if
    (
        min(materials_).value() < 0
     || max(materials_).value() > laws.size() + SMALL
    )
    {
        FatalErrorIn
        (
            "multiMaterialThermal::multiMaterialThermal\n"
            "(\n"
            "    const word& name,\n"
            "    const volScalarField& T,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Invalid definition of material indicator field.  "
            << "Number of materials: " << laws.size()
            << " max index: " << max(materials_)
            << ".  Should be " << laws.size() - 1
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiMaterialThermal::~multiMaterialThermal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::multiMaterialThermal::C() const
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
            dimensionedScalar("zeroC", dimSpecificHeatCapacity, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& result = tresult();

    // Accumulate data for all fields
    const PtrList<thermalLaw>& laws = *this;

    forAll (laws, lawI)
    {
        result.internalField() +=
            indicator(lawI)*laws[lawI].C()().internalField();
    }

    result.correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::multiMaterialThermal::k() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "k",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zerok", dimThermalConductivity, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& result = tresult();

    // Accumulate data for all fields
    const PtrList<thermalLaw>& laws = *this;

    forAll (laws, lawI)
    {
        result.internalField() +=
            indicator(lawI)*laws[lawI].k()().internalField();
    }

    result.correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::multiMaterialThermal::alpha() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "alpha",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zeroE", dimless/dimTemperature, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& result = tresult();

    // Accumulate data for all fields
    const PtrList<thermalLaw>& laws = *this;

    forAll (laws, lawI)
    {
        result.internalField() +=
            indicator(lawI)*laws[lawI].alpha()().internalField();
    }

    result.correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::multiMaterialThermal::T0() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "T0",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zeroT0", dimTemperature, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& result = tresult();

    // Accumulate data for all fields
    const PtrList<thermalLaw>& laws = *this;

    forAll (laws, lawI)
    {
        result.internalField() +=
            indicator(lawI)*laws[lawI].T0()().internalField();
    }

    result.correctBoundaryConditions();

    return tresult;
}


void Foam::multiMaterialThermal::correct()
{
    PtrList<thermalLaw>& laws = *this;

    forAll (laws, lawI)
    {
        laws[lawI].correct();
    }
}


// ************************************************************************* //
