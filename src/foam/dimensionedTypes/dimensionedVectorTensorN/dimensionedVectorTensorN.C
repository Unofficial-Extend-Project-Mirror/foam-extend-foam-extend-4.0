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

Description
    Dimensioned VectorN and TensorN obtained from generic dimensioned type.

\*---------------------------------------------------------------------------*/

#include "dimensionedVectorTensorN.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(returnType, type, fun, text)                 \
inline dimensioned< returnType > fun(const dimensioned< type >& t)  \
{                                                                   \
    return dimensioned< returnType >                                \
    (                                                               \
        #text "(" + t.name() + ')',                                 \
        fun(t.dimensions()),                                        \
        fun(t.value())                                              \
    );                                                              \
}


#define BINARY_OPERATOR(returnType, type1, type2, op, text)            \
dimensioned< returnType > op(const dimensioned< type1 >& dt1,          \
    const dimensioned< type2 >& dt2)                                   \
{                                                                      \
    return dimensioned<returnType>                                     \
    (                                                                  \
        '(' + dt1.name() + #text + dt2.name() + ')',                   \
        op(dt1.dimensions(), dt2.dimensions()),                        \
        op(dt1.value(), dt2.value())                                   \
    );                                                                 \
}

#define dimensionedType_Funs(tensorType, diagTensorType,                            \
    sphericalTensorType, vectorType, cmptType, args...)                             \
UNARY_FUNCTION(tensorType, tensorType, inv, inv)                                    \
UNARY_FUNCTION(diagTensorType, diagTensorType, inv, inv)                            \
UNARY_FUNCTION(sphericalTensorType, sphericalTensorType, inv, inv)                  \
                                                                                    \
UNARY_FUNCTION(diagTensorType, tensorType, diag, diag)                              \
UNARY_FUNCTION(diagTensorType, diagTensorType, diag, diag)                          \
UNARY_FUNCTION(sphericalTensorType, sphericalTensorType, diag, diag)                \
                                                                                    \
BINARY_OPERATOR(tensorType, tensorType, diagTensorType, operator+, +)               \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, operator+, +)               \
BINARY_OPERATOR(tensorType, tensorType, sphericalTensorType, operator+, +)          \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, operator+, +)          \
BINARY_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, operator+, +)  \
BINARY_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, operator+, +)  \
                                                                                    \
BINARY_OPERATOR(tensorType, tensorType, diagTensorType, operator-, -)               \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, operator-, -)               \
BINARY_OPERATOR(tensorType, tensorType, sphericalTensorType, operator-, -)          \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, operator-, -)          \
BINARY_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, operator-, -)  \
BINARY_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, operator-, -)

forAllVectorTensorNTypes(dimensionedType_Funs)

#undef dimensionedType_Funs
#undef UNARY_FUNCTION
#undef BINARY_OPERATOR

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
