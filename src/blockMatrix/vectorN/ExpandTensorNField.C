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
    Global functions for expansion and contraction of tensor field
    to diagonal type

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "Field.H"
#include "VectorN.H"
#include "VectorN.H"
#include "ExpandTensorN.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Cmpt, int length>
void contractScalar
(
    Field<scalar>& res,
    const UList<TensorN<Cmpt, length> >& f
)
{
    forAll (res, i)
    {
        contractScalar(res[i], f[i]);
    }
}


template <class Cmpt, int length>
void contractLinear
(
    Field<VectorN<Cmpt, length> >& res,
    const UList<TensorN<Cmpt, length> >& f
)
{
    forAll (res, i)
    {
        contractLinear(res[i], f[i]);
    }
}


template <class Cmpt, int length>
void expandScalar
(
    Field<VectorN<Cmpt, length> >& res,
    const UList<scalar>& f
)
{
    forAll (res, i)
    {
        expandScalar(res[i], f[i]);
    }
}


template <class Cmpt, int length>
void expandScalar
(
    Field<TensorN<Cmpt, length> >& res,
    const UList<scalar>& f
)
{
    forAll (res, i)
    {
        expandScalar(res[i], f[i]);
    }
}


template <class Cmpt, int length>
void expandLinear
(
    Field<TensorN<Cmpt, length> >& res,
    const UList<VectorN<Cmpt, length> >& f
)
{
    forAll (res, i)
    {
        expandLinear(res[i], f[i]);
    }
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
