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

#include "objectRegistry.H"
#include "equationReader.H"

// * * * * * * * * * * * * * scalar specializations  * * * * * * * * * * * * //

namespace Foam
{

template<>
void equationReader::evaluateTypeField
(
    Field<scalar>& resultField,
    const word& componentName,
    const word& equationName,
    const label geoIndex
) const
{
    evaluateScalarField(resultField, equationName, geoIndex);
}


template<>
void equationReader::evaluateTypeField
(
    Field<scalar>& resultField,
    const label componentIndex,
    const label equationIndex,
    const label geoIndex
) const
{
    evaluateScalarField(resultField, equationIndex, geoIndex);
}

/* DLFG 2011-09-10
 * Today I learned: Partial specialization is not allowed for member functions.
 * This code is left here to note what I want to implement in the future.
 * Apparently function objects make a nice work-around.  For now:
 *      scalar fields must use evaluate...ScalarField functions;
 *      other types must use evaluate...TypeField functions;
 * This is still safe as a mismatch produces compile-time errors.

template<class GeoMesh>
void equationReader::evaluateDimensionedTypeField<scalar, GeoMesh>
(
    DimensionedField<scalar, GeoMesh>& resultDField,
    const word& componentName,
    const word& equationName,
    const label geoIndex
) const
{
    evaluateDimensionedScalarField(resultDField, equationName, geoIndex);
}


template<class GeoMesh>
void equationReader::evaluateDimensionedTypeField<scalar, GeoMesh>
(
    DimensionedField<scalar, GeoMesh>& resultDField,
    const label componentIndex,
    const label equationIndex,
    const label geoIndex
) const
{
    evaluateDimensionedScalarField(resultDField, equationIndex, geoIndex);
}


template<template<class> class PatchField, class GeoMesh>
void equationReader::evaluateGeometricTypeField<scalar, PatchField, GeoMesh>
(
    GeometricField<scalar, PatchField, GeoMesh>& resultGField,
    const word& componentName,
    const word& equationName
) const
{
    evaluateGeometricScalarField(resultGField, equationName);
}


template<template<class> class PatchField, class GeoMesh>
void equationReader::evaluateGeometricTypeField<scalar, PatchField, GeoMesh>
(
    GeometricField<scalar, PatchField, GeoMesh>& resultGField,
    const label componentIndex,
    const label equationIndex
) const
{
    evaluateGeometricScalarField(resultGField, equationIndex);
}
*/

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
