/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Description
    Zoned multi-material rheology controlled by a material indicator field.

\*---------------------------------------------------------------------------*/

#include "multiMaterial.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiMaterial, 0);
    addToRunTimeSelectionTable(rheologyLaw, multiMaterial, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::multiMaterial::indicator
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
Foam::multiMaterial::multiMaterial
(
    const word& name,
    const volSymmTensorField& sigma,
    const dictionary& dict
)
:
    rheologyLaw(name, sigma, dict),
    PtrList<rheologyLaw>(),
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
    PtrList<rheologyLaw>& laws = *this;

    PtrList<entry> lawEntries(dict.lookup("laws"));
    laws.setSize(lawEntries.size());

    forAll (laws, lawI)
    {
        laws.set
        (
            lawI,
            rheologyLaw::New
            (
                lawEntries[lawI].keyword(),
                sigma,
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
            "multiMaterial::multiMaterial\n"
            "(\n"
            "    const word& name,\n"
            "    const volSymmTensorField& sigma,\n"
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

Foam::multiMaterial::~multiMaterial()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::multiMaterial::rho() const
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
            dimensionedScalar("zeroRho", dimDensity, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& result = tresult();

    // Accumulate data for all fields
    const PtrList<rheologyLaw>& laws = *this;

    forAll (laws, lawI)
    {
        result.internalField() +=
            indicator(lawI)*laws[lawI].rho()().internalField();
    }

    result.correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::multiMaterial::E() const
{
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
            dimensionedScalar("zeroE", dimForce/dimArea, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& result = tresult();

    // Accumulate data for all fields
    const PtrList<rheologyLaw>& laws = *this;

    forAll (laws, lawI)
    {
        result.internalField() +=
            indicator(lawI)*laws[lawI].E()().internalField();
    }

    result.correctBoundaryConditions();

    return tresult;
}


Foam::tmp<Foam::volScalarField> Foam::multiMaterial::nu() const
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
            dimensionedScalar("zeroE", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& result = tresult();

    // Accumulate data for all fields
    const PtrList<rheologyLaw>& laws = *this;

    forAll (laws, lawI)
    {
        result.internalField() +=
            indicator(lawI)*laws[lawI].nu()().internalField();
    }

    result.correctBoundaryConditions();

    return tresult;
}


void Foam::multiMaterial::correct()
{
    PtrList<rheologyLaw>& laws = *this;

    forAll (laws, lawI)
    {
        laws[lawI].correct();
    }
}


// ************************************************************************* //
