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
    Global functions for expansion and contraction of tensorN coefficient
    to diagonal and scalar type

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef expandTensorN_H
#define expandTensorN_H

#include "VectorN.H"
#include "VectorNFieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_TEMPLATE_FUNCTION(ReturnType, Type, Func)                      \
template<class Cmpt, int length>                                             \
inline void Func(ReturnType<Cmpt, length>&, const Type<Cmpt, length>&);

#define UNARY_TEMPLATE_FUNCTION_VS(ReturnType, Func)                         \
template<class Cmpt, int length>                                             \
inline void Func(ReturnType<Cmpt, length>&, const Cmpt&);

#define UNARY_TEMPLATE_FUNCTION_SV(Type, Func)                               \
template<class Cmpt, int length>                                             \
inline void Func(Cmpt&, const Type<Cmpt, length>&);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_TEMPLATE_FUNCTION_SV(TensorN, contractScalar)
UNARY_TEMPLATE_FUNCTION_SV(DiagTensorN, contractScalar)
UNARY_TEMPLATE_FUNCTION_SV(SphericalTensorN, contractScalar)
UNARY_TEMPLATE_FUNCTION_SV(VectorN, contractScalar)

UNARY_TEMPLATE_FUNCTION(VectorN, TensorN, contractLinear)
UNARY_TEMPLATE_FUNCTION(VectorN, DiagTensorN, contractLinear)
UNARY_TEMPLATE_FUNCTION(VectorN, SphericalTensorN, contractLinear)

UNARY_TEMPLATE_FUNCTION_VS(VectorN, expandScalar)
UNARY_TEMPLATE_FUNCTION_VS(TensorN, expandScalar)
UNARY_TEMPLATE_FUNCTION_VS(DiagTensorN, expandScalar)
UNARY_TEMPLATE_FUNCTION_VS(SphericalTensorN, expandScalar)

UNARY_TEMPLATE_FUNCTION(TensorN, VectorN, expandLinear)
UNARY_TEMPLATE_FUNCTION(DiagTensorN, VectorN, expandLinear)
UNARY_TEMPLATE_FUNCTION(SphericalTensorN, VectorN, expandLinear)

UNARY_TEMPLATE_FUNCTION(VectorN, TensorN, sumToDiag)
UNARY_TEMPLATE_FUNCTION(VectorN, TensorN, sumMagToDiag)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ExpandTensorNI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef UNARY_TEMPLATE_FUNCTION
#undef UNARY_TEMPLATE_FUNCTION_VS
#undef UNARY_TEMPLATE_FUNCTION_SV

#endif

// ************************************************************************* //
