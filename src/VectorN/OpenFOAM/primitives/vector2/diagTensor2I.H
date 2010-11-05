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
    diagTensor2

Description
    DiagTensorN of 2 scalars.

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Assign to a SphericalTensorN
template<>
inline void diagTensor2::operator=(const sphericalTensor2& st)
{
    this->v_[0] = st[0];
    this->v_[1] = st[0];
}

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Inner-product between two diagonal tensors
inline diagTensor2
operator&(const diagTensor2& dt1, const diagTensor2& dt2)
{
    diagTensor2 result;

    result[0] = dt1[0]*dt2[0];
    result[1] = dt1[1]*dt2[1];

    return result;
}


//- Inner-product between a spherical tensor and a diagonal tensor
inline diagTensor2
operator&(const sphericalTensor2& st1, const diagTensor2& dt2)
{
    diagTensor2 result;

    result[0] = st1[0]*dt2[0];
    result[1] = st1[0]*dt2[1];

    return result;
}

//- Inner-product between a diagonal tensor and a spherical tensor
inline diagTensor2
operator&(const diagTensor2& dt1, const sphericalTensor2& st2)
{
    diagTensor2 result;

    result[0] = dt1[0]*st2[0];
    result[1] = dt1[1]*st2[0];

    return result;
}


//- Inner-product between a diagonal tensor and a vector
inline vector2
operator&(const diagTensor2& dt, const vector2& v)
{
    vector2 result;

    result[0] = dt[0]*v[0];
    result[1] = dt[1]*v[1];

    return result;
}


//- Inner-product between a vector and a diagonal tensor
inline vector2
operator&(const vector2& v, const diagTensor2& dt)
{
    vector2 result;

    result[0] = v[0]*dt[0];
    result[1] = v[1]*dt[1];

    return result;
}


//- Division of a scalar by a diagTensor2
inline diagTensor2 operator/(const scalar s, const diagTensor2& dt)
{
    diagTensor2 result;

    result[0] = s/dt[0];
    result[1] = s/dt[1];

    return result;
}


//- Return the inverse of a diagTensor2
inline diagTensor2 inv(const diagTensor2& dt)
{
    diagTensor2 result;

    result[0] = 1.0/dt[0];
    result[1] = 1.0/dt[1];

    return result;
}


//- Return the component sum
inline scalar sum(const diagTensor2& dt)
{
    return dt[0] + dt[1];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
