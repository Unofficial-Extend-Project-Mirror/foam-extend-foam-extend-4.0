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
    GeometricTensorNFields

Description
    Specialisation of GeometricField<T> for TensorN.

SourceFiles
    GeometricTensorNFields.C

\*---------------------------------------------------------------------------*/

#ifndef GeometricTensorNFields_H
#define GeometricTensorNFields_H

#include "GeometricField.H"
#include "DimensionedTensorNFields.H"
#include "TensorNFieldFields.H"

#define TEMPLATE template<template<class> class PatchField, class GeoMesh>
#include "GeometricFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define TensorN_FieldFunctions(tensorType, diagTensorType, sphericalTensorType,     \
    vectorType, CmptType, args...)                                                  \
                                                                                    \
UNARY_FUNCTION(tensorType, tensorType,inv,inv)                                      \
UNARY_FUNCTION(diagTensorType, tensorType,diag,diag)                                \
UNARY_FUNCTION(tensorType, tensorType,negSumDiag,negSumDiag)                        \
UNARY_FUNCTION(CmptType, tensorType,contractScalar,contractScalar)                  \
UNARY_FUNCTION(vectorType, tensorType,contractLinear,contractLinear)                \
                                                                                    \
BINARY_OPERATOR(tensorType, CmptType, tensorType, *,'*',multiply)                   \
BINARY_TYPE_OPERATOR(tensorType, CmptType, tensorType, *,'*',multiply)              \
                                                                                    \
BINARY_OPERATOR(tensorType, CmptType, tensorType, /,'|',divide)                     \
BINARY_TYPE_OPERATOR(tensorType, CmptType, tensorType, /,'|',divide)                \
                                                                                    \
BINARY_OPERATOR(vectorType, vectorType, tensorType, /,'|',divide)                   \
BINARY_TYPE_OPERATOR(vectorType, vectorType, tensorType, /,'|',divide)              \
                                                                                    \
BINARY_OPERATOR(tensorType, tensorType, tensorType, /,'|',divide)                   \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, /,'|',divide)              \
                                                                                    \
BINARY_OPERATOR(tensorType, tensorType, diagTensorType, /,'|',divide)               \
BINARY_TYPE_OPERATOR(tensorType, tensorType, diagTensorType, /,'|',divide)          \
                                                                                    \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, /,'|',divide)               \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, /,'|',divide)          \
                                                                                    \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, /,'|',divide)          \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, /,'|',divide)     \
                                                                                    \
BINARY_OPERATOR(tensorType, tensorType, sphericalTensorType, /,'|',divide)          \
BINARY_TYPE_OPERATOR(tensorType, tensorType, sphericalTensorType, /,'|',divide)     \
                                                                                    \
BINARY_OPERATOR(tensorType, tensorType, tensorType, +,'+',add)                      \
BINARY_OPERATOR(tensorType, tensorType, tensorType, -,'-',subtract)                 \
                                                                                    \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, +,'+',add)                 \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, -,'-',subtract)            \
                                                                                    \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, +,'+',add)                  \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, -,'-',subtract)             \
                                                                                    \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, +,'+',add)             \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, -,'-',subtract)        \
                                                                                    \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, +,'+',add)             \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, -,'-',subtract)        \
                                                                                    \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, +,'+',add)        \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, -,'-',subtract)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forAllVectorTensorNTypes(TensorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef TensorN_FieldFunctions

#include "undefFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "GeometricTensorNFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
