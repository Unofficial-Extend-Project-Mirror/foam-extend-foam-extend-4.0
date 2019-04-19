/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "multiMaterialThermal.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiMaterialThermal, 0);
    addToRunTimeSelectionTable(thermalLaw, multiMaterialThermal, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::multiMaterialThermal::indicator
(
    const label i
) const
{
    const scalarField& mat = materials_.internalField();

    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "indicator",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimless,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& result = tresult();

    forAll (mat, matI)
    {
        if (mat[matI] > i - SMALL && mat[matI] < i + SMALL)
        {
            result[matI] = 1.0;
        }
        else
        {
            result[matI] = 0.0;
        }
    }

    result.correctBoundaryConditions();

    return tresult;
}


void Foam::multiMaterialThermal::readLaws
(
    const volScalarField& T,
    const dictionary& dict
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
}


void Foam::multiMaterialThermal::checkLaws() const
{
    const PtrList<thermalLaw>& laws = *this;

    if
    (
        max(materials_).value() > laws.size() + SMALL
    )
    {
        FatalErrorIn
        (
            "multiMaterialThermal::checkLaws()\n"
        )   << "Invalid definition of material indicator field.  "
            << "Number of materials: " << laws.size()
            << " max index: " << max(materials_)
            << ".  Should be " << laws.size() - 1
            << abort(FatalError);
    }

    if
    (
        min(materials_).value() < 0
    )
    {
        FatalErrorIn
        (
            "multiMaterialThermal::checkLaws()\n"
        )   << "Invalid definition of material indicator field.  "
            << "Number of materials: " << laws.size()
            << " min index: " << min(materials_)
            << ".  Should be 0"
            << abort(FatalError);
    }
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
    readLaws(T, dict);
    checkLaws();
}


// Construct from dictionary and create default material field
Foam::multiMaterialThermal::multiMaterialThermal
(
    const word& name,
    const volScalarField& T,
    const dictionary& dict,
    const scalar defaultMaterial
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
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        defaultMaterial,
        zeroGradientFvPatchScalarField::typeName
    )
{
    readLaws(T, dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::multiMaterialThermal::~multiMaterialThermal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::multiMaterialThermal::rho() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "rhoTmp",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zeroRho", dimDensity, 0),
            calculatedFvPatchScalarField::typeName
        )
    );
    volScalarField& result = tresult();

    // Accumulate data for all fields
    const PtrList<thermalLaw>& laws = *this;

    forAll (laws, lawI)
    {
        result += indicator(lawI)*laws[lawI].rho()();
    }

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::multiMaterialThermal::C() const
{
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
                "CTmp",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zeroC", dimSpecificHeatCapacity, 0),
            calculatedFvPatchScalarField::typeName
        )
    );
    volScalarField& result = tresult();

    // Accumulate data for all fields
    const PtrList<thermalLaw>& laws = *this;

    forAll (laws, lawI)
    {
        result += indicator(lawI)*laws[lawI].C()();
    }

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
                "kTmp",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zerok", dimensionSet(1, 1, -3, -1, 0), 0),
            calculatedFvPatchScalarField::typeName
        )
    );
    volScalarField& result = tresult();

    // Accumulate data for all fields
    const PtrList<thermalLaw>& laws = *this;

    forAll (laws, lawI)
    {
        result += indicator(lawI)*laws[lawI].k();
    }

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
