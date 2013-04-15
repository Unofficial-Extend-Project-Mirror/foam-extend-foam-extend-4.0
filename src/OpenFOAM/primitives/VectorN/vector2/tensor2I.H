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
    tensor2

Description
    TensorN of 2 scalars.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return tensor transpose
template<>
inline tensor2 tensor2::T() const
{
    tensor2 transpose;

    transpose[0] = this->operator[](0);
    transpose[1] = this->operator[](2);
    transpose[2] = this->operator[](1);
    transpose[3] = this->operator[](3);

    return transpose;
}


//- Assign to a sphericalTensor2
template<>
inline void tensor2::operator=(const sphericalTensor2& st)
{
    this->v_[0] = st[0];
    this->v_[1] = 0.0;
    this->v_[2] = 0.0;
    this->v_[3] = st[0];
}


//- Assign to a diagTensor2
template<>
inline void tensor2::operator=(const diagTensor2& dt)
{
    this->v_[0] = dt[0];
    this->v_[1] = 0.0;
    this->v_[2] = 0.0;
    this->v_[3] = dt[1];
}

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Inner-product between two tensors
inline tensor2
operator&(const tensor2& t1, const tensor2& t2)
{
    tensor2 result;

    result[0] = t1[0]*t2[0] + t1[1]*t2[2];
    result[1] = t1[0]*t2[1] + t1[1]*t2[3];
    result[2] = t1[2]*t2[0] + t1[3]*t2[2];
    result[3] = t1[2]*t2[1] + t1[3]*t2[3];

    return result;
}


//- Inner-product between a diagonal tensors and a tensor
inline tensor2
operator&(const diagTensor2& dt1, const tensor2& t2)
{
    tensor2 result;

    result[0] = dt1[0]*t2[0];
    result[1] = dt1[0]*t2[1];
    result[2] = dt1[1]*t2[2];
    result[3] = dt1[1]*t2[3];

    return result;
}

//- Inner-product between a tensor and diagonal tensor
inline tensor2
operator&(const tensor2& t1, const diagTensor2& dt2)
{
    tensor2 result;

    result[0] = t1[0]*dt2[0];
    result[1] = t1[1]*dt2[1];
    result[2] = t1[2]*dt2[0];
    result[3] = t1[3]*dt2[1];

    return result;
}


//- Inner-product between a spherical tensor and a tensor
inline tensor2
operator&(const sphericalTensor2& st1, const tensor2& t2)
{
    tensor2 result;

    result[0] = st1[0]*t2[0];
    result[1] = st1[0]*t2[1];
    result[2] = st1[0]*t2[2];
    result[3] = st1[0]*t2[3];

    return result;
}

//- Inner-product between a tensor and spherical tensor
inline tensor2
operator&(const tensor2& t1, const sphericalTensor2& st2)
{
    tensor2 result;

    result[0] = t1[0]*st2[0];
    result[1] = t1[1]*st2[0];
    result[2] = t1[2]*st2[0];
    result[3] = t1[3]*st2[0];

    return result;
}


//- Inner-product between a tensor and a vector
inline vector2
operator&(const tensor2& t, const vector2& v)
{
    vector2 result;

    result[0] = t[0]*v[0] + t[1]*v[1];
    result[1] = t[2]*v[0] + t[3]*v[1];

    return result;
}


//- Inner-product between a vector and a tensor
inline vector2
operator&(const vector2& v, const tensor2& t)
{
    vector2 result;

    result[0] = v[0]*t[0] + v[1]*t[2];
    result[1] = v[0]*t[1] + v[1]*t[3];

    return result;
}


//- Outer-product between two vectors
inline tensor2
operator*(const vector2& v1, const vector2& v2)
{
    tensor2 result;

    result[0] = v1[0]*v2[0];
    result[1] = v1[0]*v2[1];
    result[2] = v1[1]*v2[0];
    result[3] = v1[1]*v2[1];

    return result;
}


//- Return the determinant of a tensor
inline scalar det(const tensor2& t)
{
    return
    (
        t[0]*t[3]-t[1]*t[2]
    );
}

//- Return the inverse of a tensor given the determinant
inline tensor2 inv(const tensor2& t)
{
    tensor2 cofactor;

    cofactor[0] = t[3];
    cofactor[1] = -t[1];
    cofactor[2] = -t[2];
    cofactor[3] = t[0];

    return cofactor/det(t);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
