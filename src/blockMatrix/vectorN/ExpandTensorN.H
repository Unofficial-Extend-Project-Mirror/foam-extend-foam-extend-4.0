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

Description
    Global functions for expansion and contraction of tensor8 coefficient
    to diagonal type

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef expandTensorN_H
#define expandTensorN_H

#include "VectorN.H"
#include "TensorN.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Return the diagonal of a TensorN as a scalar
template <class Cmpt, int length>
inline void contractScalar(Cmpt& result, const TensorN<Cmpt, length>& t)
{
    result = pTraits<Cmpt>::zero;

    for (int i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        result += t(i, i);
    }

    result *= 1.0/8.0;
}


//- Return the diagonal of a TensorN as a vector
template <class Cmpt, int length>
inline void contractLinear
(
    VectorN<Cmpt, length>& result,
    const TensorN<Cmpt, length>& t
)
{
    result = VectorN<Cmpt, length>::zero;

    for (int i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        result(i) = t(i, i);
    }
}


//- Return the VectorN given a scalar
template <class Cmpt, int length>
inline void expandScalar(VectorN<Cmpt, length>& result, const Cmpt& v)
{
    result = VectorN<Cmpt, length>::zero;

    for (int i = 0; i < VectorN<Cmpt, length>::nComponents; i++)
    {
        result(i) = v;
    }
}


//- Return the TensorN given a scalar
template <class Cmpt, int length>
inline void expandScalar(TensorN<Cmpt, length>& result, const Cmpt& v)
{
    result = TensorN<Cmpt, length>::zero;

    for (int i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        result(i, i) = v;
    }
}


//- Return the TensorN given a diagonal vector8
template <class Cmpt, int length>
inline void expandLinear
(
    TensorN<Cmpt, length>& result,
    const VectorN<Cmpt, length>& v
)
{
    result = TensorN<Cmpt, length>::zero;

    for (int i = 0; i < TensorN<Cmpt, length>::rowLength; i++)
    {
        result(i, i) = v(i);
    }
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
