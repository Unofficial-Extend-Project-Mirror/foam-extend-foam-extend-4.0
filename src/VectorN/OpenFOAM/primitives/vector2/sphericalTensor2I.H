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

Type
    sphericalTensor2

Description
    SphericalTensorN of 2 scalars.

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Inner-product between two spherical tensors
inline sphericalTensor2
operator&(const sphericalTensor2& st1, const sphericalTensor2& st2)
{
    return sphericalTensor2(st1[0]*st2[0]);
}


//- Inner-product between a spherical tensor and a vector
inline vector2
operator&(const sphericalTensor2& st, const vector2& v)
{
    vector2 result;

    result[0] = st[0]*v[0];
    result[1] = st[0]*v[1];

    return result;
}


//- Inner-product between a vector and a spherical tensor
inline vector2
operator&(const vector2& v, const sphericalTensor2& st)
{
    vector2 result;

    result[0] = v[0]*st[0];
    result[1] = v[1]*st[0];

    return result;
}


//- Return the inverse of a tensor given the determinant
inline sphericalTensor2 inv(const sphericalTensor2& st)
{
    return sphericalTensor2(1/st[0]);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
