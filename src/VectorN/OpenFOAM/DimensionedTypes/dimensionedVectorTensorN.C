/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

Description
    Dimensioned VectorN and TensorN obtained from generic dimensioned type.

\*---------------------------------------------------------------------------*/

#include "dimensionedVectorTensorN.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(returnType, type, fun, text)                 \
inline returnType fun(const type& t)                                \
{                                                                   \
    return returnType                                               \
    (                                                               \
        #text "(" + t.name() + ')',                                 \
        fun(t.dimensions()),                                        \
        fun(t.value())                                              \
    );                                                              \
}

#define dimensionedType_Funs(cmptType, vectorType, tensorType,      \
    diagTensorType, sphericalTensorType)                            \
UNARY_FUNCTION(tensorType, tensorType, inv, inv)                    \
UNARY_FUNCTION(diagTensorType, diagTensorType, inv, inv)            \
UNARY_FUNCTION(sphericalTensorType, sphericalTensorType, inv, inv)  \
UNARY_FUNCTION(diagTensorType, tensorType, diag, diag)              \
UNARY_FUNCTION(diagTensorType, diagTensorType, diag, diag)          \
UNARY_FUNCTION(sphericalTensorType, sphericalTensorType, diag, diag)\
UNARY_FUNCTION(cmptType, vectorType, cmptSum, cmptSum)

#define dimensionedVectorN_Funs(length)                         \
dimensionedType_Funs                                            \
(                                                               \
    dimensionedScalar,                                          \
    dimensionedVector##length,                                  \
    dimensionedTensor##length,                                  \
    dimensionedDiagTensor##length,                              \
    dimensionedSphericalTensor##length                          \
)

dimensionedVectorN_Funs(2)
dimensionedVectorN_Funs(4)
dimensionedVectorN_Funs(6)
dimensionedVectorN_Funs(8)

#undef dimensionedVectorN_Funs
#undef dimensionedType_Funs
#undef UNARY_FUNCTION

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
