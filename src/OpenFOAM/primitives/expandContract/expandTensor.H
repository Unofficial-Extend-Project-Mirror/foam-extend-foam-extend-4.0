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
    Global functions for expansion and contraction of tensor coefficient
    to diagonal type

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef expandTensor_H
#define expandTensor_H

#include "Tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Return the diagonal of a tensor as a scalar
template <class Cmpt>
inline void contractScalar(Cmpt& result, const Tensor<Cmpt>& t)
{
    result = 1.0/3.0*(t.xx() + t.yy() + t.zz());
}


//- Return the diagonal of a tensor as a scalar
template <class Cmpt>
inline Cmpt contractScalar(const Tensor<Cmpt>& t)
{
    Cmpt result;
    contractScalar(result, t);
    return result;
}


//- Return the diagonal of a tensor as a vector
template <class Cmpt>
inline void contractLinear(Vector<Cmpt>& result, const Tensor<Cmpt>& t)
{
    result = Vector<Cmpt>(t.xx(), t.yy(), t.zz());
}


//- Return the diagonal of a tensor as a vector
template <class Cmpt>
inline Vector<Cmpt> contractLinear(const Tensor<Cmpt>& t)
{
    Vector<Cmpt> result;
    contractLinear(result, t);
    return result;
}


//- Return the vector given a scalar
template <class Cmpt>
inline void expandScalar(Vector<Cmpt>& result, const Cmpt& v)
{
    result = Vector<Cmpt>(v, v, v);
}


//- Return the tensor given a scalar
template <class Cmpt>
inline void expandScalar(Tensor<Cmpt>& result, const Cmpt& v)
{
    result = Tensor<Cmpt>
    (
        v,                    pTraits<Cmpt>::zero,     pTraits<Cmpt>::zero,
        pTraits<Cmpt>::zero,  v,                       pTraits<Cmpt>::zero,
        pTraits<Cmpt>::zero,  pTraits<Cmpt>::zero,     v
    );
}


//- Return the tensor given a diagonal vector
template <class Cmpt>
inline void expandLinear(Tensor<Cmpt>& result, const Vector<Cmpt>& v)
{
    result = Tensor<Cmpt>
    (
        v.x(),                pTraits<Cmpt>::zero,     pTraits<Cmpt>::zero,
        pTraits<Cmpt>::zero,  v.y(),                   pTraits<Cmpt>::zero,
        pTraits<Cmpt>::zero,  pTraits<Cmpt>::zero,     v.z()
    );
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
