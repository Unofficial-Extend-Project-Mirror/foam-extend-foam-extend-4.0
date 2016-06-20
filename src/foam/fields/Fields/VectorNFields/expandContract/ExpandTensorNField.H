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
    Global functions for expansion and contraction of tensor field
    to diagonal type

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    ExpandTensorNField.C

\*---------------------------------------------------------------------------*/

#ifndef expandTensorNField_H
#define expandTensorNField_H

#include "Field.H"
#include "VectorNFieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(typeF1, typeF2, FUNC)                                  \
void FUNC(Field<typeF1>& f1, const UList<typeF2>& f2);                        \
void FUNC(Field<typeF1>& f1, const tmp<Field<typeF2> >& tf2);


#define ExpandFieldFunctions(tensorType, diagTensorType, sphericalTensorType, \
        vectorType, cmptType, args...)                                        \
                                                                              \
UNARY_FUNCTION(cmptType, tensorType, contractScalar)                          \
UNARY_FUNCTION(cmptType, diagTensorType, contractScalar)                      \
UNARY_FUNCTION(cmptType, sphericalTensorType, contractScalar)                 \
UNARY_FUNCTION(cmptType, vectorType, contractScalar)                          \
                                                                              \
UNARY_FUNCTION(vectorType, tensorType, contractLinear)                        \
UNARY_FUNCTION(vectorType, diagTensorType, contractLinear)                    \
UNARY_FUNCTION(vectorType, sphericalTensorType, contractLinear)               \
                                                                              \
UNARY_FUNCTION(vectorType, cmptType, expandScalar)                            \
UNARY_FUNCTION(tensorType, cmptType, expandScalar)                            \
UNARY_FUNCTION(diagTensorType, cmptType, expandScalar)                        \
UNARY_FUNCTION(sphericalTensorType, cmptType, expandScalar)                   \
                                                                              \
UNARY_FUNCTION(tensorType, vectorType, expandLinear)                          \
UNARY_FUNCTION(diagTensorType, vectorType, expandLinear)                      \
UNARY_FUNCTION(sphericalTensorType, vectorType, expandLinear)                 \
                                                                              \
UNARY_FUNCTION(vectorType, tensorType, sumToDiag)                             \
UNARY_FUNCTION(vectorType, tensorType, sumMagToDiag)


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

namespace Foam
{

forAllVectorTensorNTypes(ExpandFieldFunctions)

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef UNARY_FUNCTION
#undef ExpandFieldFunctions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "ExpandTensorNField.C"
// #endif

// ************************************************************************* //

#endif
