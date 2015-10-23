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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::equationReader::evaluateTypeField
(
    Field<Type>& resultField,
    const word& componentName,
    const word& equationName,
    const label geoIndex
) const
{
    // Get equationIndex
    label equationIndex(lookup(equationName));
    if (equationIndex < 0)
    {
        FatalErrorIn("equationReader::evaluateTypeField")
            << "Equation name " << equationName << " not found."
            << abort(FatalError);
    }

    // Get componentIndex
    Type dummy;
    label componentIndex(-1);
    forAll(dummy, i)
    {
        if (Type::componentNames[i] == componentName)
        {
            componentIndex = i;
            break;
        }
    }
    if (componentIndex < 0)
    {
        wordList validNames(dummy.size());
        forAll(dummy, i)
        {
            validNames[i] = Type::componentNames[i];
        }
        FatalErrorIn("equationReader::evaluateTypeField")
            << componentName << " is not a valid component name.  Valid names "
            << "are " << validNames
            << abort(FatalError);
    }

    evaluateTypeField(resultField, componentIndex, equationIndex, geoIndex);
}


template<class Type>
void Foam::equationReader::evaluateTypeField
(
    Field<Type>& resultField,
    const label componentIndex,
    const label equationIndex,
    const label geoIndex
) const
{
    scalarField sf(resultField.size());
    evaluateScalarField(sf, equationIndex, geoIndex);

    List_ACCESS(Type, resultField, resultFieldP);
    List_CONST_ACCESS(scalar, sf, sfP);

    /* loop through fields performing f1 OP1 f2 OP2 f3 */
    List_FOR_ALL(resultField, i)
        List_ELEM(resultField, resultFieldP, i)[componentIndex] =
            List_ELEM(sf, sfP, i);
    List_END_FOR_ALL
}


template<class GeoMesh>
void Foam::equationReader::evaluateDimensionedScalarField
(
    DimensionedField<scalar, GeoMesh>& resultDField,
    const word& equationName,
    const label geoIndex
) const
{
    label equationIndex(lookup(equationName));
    if (equationIndex < 0)
    {
        FatalErrorIn("equationReader::evaluateDimensionedScalarField")
            << "Equation name " << equationName << " not found."
            << abort(FatalError);
    }
    evaluateDimensionedScalarField(resultDField, equationIndex, geoIndex);
}


template<class GeoMesh>
void Foam::equationReader::evaluateDimensionedScalarField
(
    DimensionedField<scalar, GeoMesh>& resultDField,
    const label equationIndex,
    const label geoIndex
) const
{
    evaluateScalarField(resultDField.field(), equationIndex, geoIndex);
    checkFinalDimensions
    (
        equationIndex,
        resultDField.dimensions(),
        resultDField.name()
    );
}


template<class Type, class GeoMesh>
void Foam::equationReader::evaluateDimensionedTypeField
(
    DimensionedField<Type, GeoMesh>& resultDField,
    const word& componentName,
    const word& equationName,
    const label geoIndex
) const
{
    // Get equationIndex
    label equationIndex(lookup(equationName));
    if (equationIndex < 0)
    {
        FatalErrorIn("equationReader::evaluateDimensionedTypeField")
            << "Equation name " << equationName << " not found."
            << abort(FatalError);
    }

    // Get componentIndex
    Type dummy;
    label componentIndex(-1);
    forAll(dummy, i)
    {
        if (Type::componentNames[i] == componentName)
        {
            componentIndex = i;
            break;
        }
    }
    if (componentIndex < 0)
    {
        wordList validNames(dummy.size());
        forAll(dummy, i)
        {
            validNames[i] = Type::componentNames[i];
        }
        FatalErrorIn("equationReader::evaluateDimensionedTypeField")
            << componentName << " is not a valid component name.  Valid names "
            << "are " << validNames
            << abort(FatalError);
    }

    evaluateDimensionedTypeField
    (
        resultDField,
        componentIndex,
        equationIndex,
        geoIndex
    );
}


template<class Type, class GeoMesh>
void Foam::equationReader::evaluateDimensionedTypeField
(
    DimensionedField<Type, GeoMesh>& resultDField,
    const label componentIndex,
    const label equationIndex,
    const label geoIndex
) const
{
    evaluateTypeField
    (
        resultDField.field(),
        componentIndex,
        equationIndex,
        geoIndex
    );
    checkFinalDimensions
    (
        equationIndex,
        resultDField.dimensions(),
        resultDField.name() + "." + Type::componentNames[componentIndex]
    );
}


template<template<class> class PatchField, class GeoMesh>
void Foam::equationReader::evaluateGeometricScalarField
(
    GeometricField<scalar, PatchField, GeoMesh>& resultGField,
    const word& equationName
) const
{
    label equationIndex(lookup(equationName));
    if (equationIndex < 0)
    {
        FatalErrorIn("equationReader::evaluateGeometricScalarField")
            << "Equation name " << equationName << " not found."
            << abort(FatalError);
    }
    evaluateGeometricScalarField(resultGField, equationIndex);
}


template <template<class> class PatchField, class GeoMesh>
void Foam::equationReader::evaluateGeometricScalarField
(
    GeometricField<scalar, PatchField, GeoMesh>& resultGField,
    const label equationIndex
) const
{
    // Internal field: geoIndex = 0
    evaluateScalarField
    (
        resultGField.internalField(),
        equationIndex,
        0
    );
    // Boundary fields: geoIndex = patchIndex + 1
    forAll(resultGField.boundaryField(), patchIndex)
    {
        evaluateScalarField
        (
            resultGField.boundaryField()[patchIndex],
            equationIndex,
            patchIndex + 1
        );
    }
    checkFinalDimensions
    (
        equationIndex,
        resultGField.dimensions(),
        resultGField.name()
    );
}


template
<
    class Type, template<class> class PatchField, class GeoMesh
>
void Foam::equationReader::evaluateGeometricTypeField
(
    GeometricField<Type, PatchField, GeoMesh>& resultGField,
    const word& componentName,
    const word& equationName
) const
{
    // Get equationIndex
    label equationIndex(lookup(equationName));
    if (equationIndex < 0)
    {
        FatalErrorIn("equationReader::evaluateGeometricTypeField")
            << "Equation name " << equationName << " not found."
            << abort(FatalError);
    }

    // Get componentIndex
    Type dummy;
    label componentIndex(-1);
    forAll(dummy, i)
    {
        if (Type::componentNames[i] == componentName)
        {
            componentIndex = i;
            break;
        }
    }
    if (componentIndex < 0)
    {
        wordList validNames(dummy.size());
        forAll(dummy, i)
        {
            validNames[i] = Type::componentNames[i];
        }
        FatalErrorIn("equationReader::evaluateGeometricTypeField")
            << componentName << " is not a valid component name.  Valid names "
            << "are " << validNames
            << abort(FatalError);
    }
    evaluateGeometricTypeField(resultGField, componentIndex, equationIndex);
}


template
<
    class Type, template<class> class PatchField, class GeoMesh
>
void Foam::equationReader::evaluateGeometricTypeField
(
    GeometricField<Type, PatchField, GeoMesh>& resultGField,
    const label componentIndex,
    const label equationIndex
) const
{
    // Internal field: geoIndex = 0
    evaluateTypeField
    (
        resultGField.internalField(),
        componentIndex,
        equationIndex,
        0
    );

    // Boundary fields: geoIndex = patchIndex + 1
    forAll(resultGField.boundaryField(), patchIndex)
    {
        evaluateTypeField
        (
            resultGField.boundaryField()[patchIndex],
            componentIndex,
            equationIndex,
            patchIndex + 1
        );
    }
    checkFinalDimensions
    (
        equationIndex,
        resultGField.dimensions(),
        resultGField.name() + "." + Type::componentNames[componentIndex]
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
