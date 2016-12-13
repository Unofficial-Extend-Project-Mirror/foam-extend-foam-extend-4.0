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

Class
    tensorNFieldFields

Description
    Specialisation of FieldField<Field,T> for TensorN.

SourceFiles
    TensorNFieldFields.C

\*---------------------------------------------------------------------------*/

#ifndef TensorNFieldFields_H
#define TensorNFieldFields_H

#include "VectorTensorNFieldFieldsFwd.H"
#include "FieldField.H"

#define TEMPLATE template<template<class> class Field>
#include "FieldFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define TensorN_FieldFunctions(tensorType, diagTensorType, sphericalTensorType,         \
    vectorType, CmptType, args...)                                                      \
                                                                                        \
UNARY_FUNCTION(tensorType, tensorType, inv)                                             \
UNARY_FUNCTION(diagTensorType, tensorType, diag)                                        \
UNARY_FUNCTION(tensorType, tensorType, negSumDiag)                                      \
UNARY_FUNCTION(vectorType, tensorType, contractLinear)                                  \
UNARY_FUNCTION(CmptType, tensorType, contractScalar)                                    \
                                                                                        \
BINARY_OPERATOR(tensorType, CmptType, tensorType, /, divide)                            \
BINARY_TYPE_OPERATOR(tensorType, CmptType, tensorType, /, divide)                       \
                                                                                        \
BINARY_OPERATOR(vectorType, vectorType, tensorType, /, divide)                          \
BINARY_TYPE_OPERATOR(vectorType, vectorType, tensorType, /, divide)                     \
                                                                                        \
BINARY_OPERATOR(tensorType, tensorType, tensorType, /, divide)                          \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, /, divide)                     \
                                                                                        \
BINARY_OPERATOR(tensorType, tensorType, diagTensorType, /, divide)                      \
BINARY_TYPE_OPERATOR(tensorType, tensorType, diagTensorType, /, divide)                 \
                                                                                        \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, /, divide)                      \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, /, divide)                 \
                                                                                        \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, /, divide)                 \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, /, divide)            \
                                                                                        \
BINARY_OPERATOR(tensorType, tensorType, sphericalTensorType, /, divide)                 \
BINARY_TYPE_OPERATOR(tensorType, tensorType, sphericalTensorType, /, divide)            \
                                                                                        \
BINARY_OPERATOR(tensorType, tensorType, tensorType, +, add)                             \
BINARY_OPERATOR(tensorType, tensorType, tensorType, -, subtract)                        \
                                                                                        \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, +, add)                        \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, -, subtract)                   \
                                                                                        \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, +, add)                         \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, -, subtract)                    \
                                                                                        \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, +, add)                    \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, -, subtract)               \
                                                                                        \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, +, add)                    \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, -, subtract)               \
                                                                                        \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, +, add)               \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, -, subtract)


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
#   include "TensorNFieldFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
