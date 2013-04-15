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

Class
    DimensionedVectorNFields

Description
    Specialisations of DimensionedField<T> for VectorN types.

SourceFiles
    DimensionedVectorNFields.C

\*---------------------------------------------------------------------------*/

#ifndef DimensionedVectorNFields_H
#define DimensionedVectorNFields_H

#include "DimensionedField.H"
#include "dimensionedVectorTensorN.H"
#include "VectorTensorNFields.H"

#define TEMPLATE template<class GeoMesh>
#include "DimensionedFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define VectorN_FieldFunctions(tensorType, diagTensorType, sphericalTensorType, \
    vectorType, CmptType, args...)                                              \
                                                                                \
UNARY_FUNCTION(CmptType, vectorType, cmptSum, cmptSum)                          \
                                                                                \
BINARY_FUNCTION(vectorType, vectorType, vectorType, cmptMultiply)               \
BINARY_TYPE_FUNCTION(vectorType, vectorType, vectorType, cmptMultiply)          \
                                                                                \
BINARY_OPERATOR(vectorType, CmptType, vectorType, /,'|',divide)                 \
BINARY_TYPE_OPERATOR(vectorType, CmptType, vectorType, /,'|',divide)            \
                                                                                \
BINARY_OPERATOR(vectorType, vectorType, vectorType, +,'+',add)                  \
BINARY_OPERATOR(vectorType, vectorType, vectorType, -,'-',subtract)             \
                                                                                \
BINARY_TYPE_OPERATOR(vectorType, vectorType, vectorType, +,'+',add)             \
BINARY_TYPE_OPERATOR(vectorType, vectorType, vectorType, -,'-',subtract)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forAllVectorTensorNTypes(VectorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef VectorN_FieldFunctions

#include "undefFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DimensionedVectorNFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
