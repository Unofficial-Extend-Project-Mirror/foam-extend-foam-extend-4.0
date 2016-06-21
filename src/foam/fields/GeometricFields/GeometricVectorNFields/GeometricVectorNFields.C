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

#include "GeometricVectorNFields.H"

#define TEMPLATE template<template<class> class PatchField, class GeoMesh>
#include "GeometricFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define VectorN_FieldFunctions(tensorType, diagTensorType, sphericalTensorType, \
    vectorType, CmptType, args...)                                              \
                                                                                \
UNARY_FUNCTION(CmptType, vectorType, cmptSum, cmptSum)                          \
                                                                                \
BINARY_FUNCTION(vectorType, vectorType, vectorType, cmptMultiply)               \
BINARY_TYPE_FUNCTION(vectorType, vectorType, vectorType, cmptMultiply)          \
                                                                                \
BINARY_OPERATOR(vectorType, CmptType, vectorType, *,'*',multiply)               \
BINARY_TYPE_OPERATOR(vectorType, CmptType, vectorType, *,'*',multiply)          \
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

// ************************************************************************* //
