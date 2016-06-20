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

Type
    tensor4

Description
    TensorN of 4 scalars.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


//- Return the determinant of a tensor
inline scalar det(const tensor4& t)
{
    return
    (
        t[3]*t[6]*t[9]*t[12]  - t[2]*t[7]*t[9]*t[12]  - t[3]*t[5]*t[10]*t[12] + t[1]*t[7]*t[10]*t[12] +
        t[2]*t[5]*t[11]*t[12] - t[1]*t[6]*t[11]*t[12] - t[3]*t[6]*t[8]*t[13]  + t[2]*t[7]*t[8]*t[13]  +
        t[3]*t[4]*t[10]*t[13] - t[0]*t[7]*t[10]*t[13] - t[2]*t[4]*t[11]*t[13] + t[0]*t[6]*t[11]*t[13] +
        t[3]*t[5]*t[8]*t[14]  - t[1]*t[7]*t[8]*t[14]  - t[3]*t[4]*t[9]*t[14]  + t[0]*t[7]*t[9]*t[14]  +
        t[1]*t[4]*t[11]*t[14] - t[0]*t[5]*t[11]*t[14] - t[2]*t[5]*t[8]*t[15]  + t[1]*t[6]*t[8]*t[15]  +
        t[2]*t[4]*t[9]*t[15]  - t[0]*t[6]*t[9]*t[15]  - t[1]*t[4]*t[10]*t[15] + t[0]*t[5]*t[10]*t[15]
    );
}


//- Return the inverse of a tensor given the determinant
inline tensor4 inv(const tensor4& t)
{
    tensor4 cofactor;

    cofactor[0] =   t[6]*t[11]*t[13] - t[7]*t[10]*t[13] + t[7]*t[9]*t[14]
                  - t[5]*t[11]*t[14] - t[6]*t[9]*t[15]  + t[5]*t[10]*t[15];
    cofactor[1] =   t[3]*t[10]*t[13] - t[2]*t[11]*t[13] - t[3]*t[9]*t[14]
                  + t[1]*t[11]*t[14] + t[2]*t[9]*t[15]  - t[1]*t[10]*t[15];
    cofactor[2] =   t[2]*t[7]*t[13]  - t[3]*t[6]*t[13]  + t[3]*t[5]*t[14]
                  - t[1]*t[7]*t[14]  - t[2]*t[5]*t[15]  + t[1]*t[6]*t[15];
    cofactor[3] =   t[3]*t[6]*t[9]   - t[2]*t[7]*t[9]   - t[3]*t[5]*t[10]
                  + t[1]*t[7]*t[10]  + t[2]*t[5]*t[11]  - t[1]*t[6]*t[11];
    cofactor[4] =   t[7]*t[10]*t[12] - t[6]*t[11]*t[12] - t[7]*t[8]*t[14]
                  + t[4]*t[11]*t[14] + t[6]*t[8]*t[15]  - t[4]*t[10]*t[15];
    cofactor[5] =   t[2]*t[11]*t[12] - t[3]*t[10]*t[12] + t[3]*t[8]*t[14]
                  - t[0]*t[11]*t[14] - t[2]*t[8]*t[15]  + t[0]*t[10]*t[15];
    cofactor[6] =   t[3]*t[6]*t[12]  - t[2]*t[7]*t[12]  - t[3]*t[4]*t[14]
                  + t[0]*t[7]*t[14]  + t[2]*t[4]*t[15]  - t[0]*t[6]*t[15];
    cofactor[7] =   t[2]*t[7]*t[8]   - t[3]*t[6]*t[8]   + t[3]*t[4]*t[10]
                  - t[0]*t[7]*t[10]  - t[2]*t[4]*t[11]  + t[0]*t[6]*t[11];
    cofactor[8] =   t[5]*t[11]*t[12] - t[7]*t[9]*t[12]  + t[7]*t[8]*t[13]
                  - t[4]*t[11]*t[13] - t[5]*t[8]*t[15]  + t[4]*t[9]*t[15];
    cofactor[9] =   t[3]*t[9]*t[12]  - t[1]*t[11]*t[12] - t[3]*t[8]*t[13]
                  + t[0]*t[11]*t[13] + t[1]*t[8]*t[15]  - t[0]*t[9]*t[15];
    cofactor[10] =  t[1]*t[7]*t[12]  - t[3]*t[5]*t[12]  + t[3]*t[4]*t[13]
                  - t[0]*t[7]*t[13]  - t[1]*t[4]*t[15]  + t[0]*t[5]*t[15];
    cofactor[11] =  t[3]*t[5]*t[8]   - t[1]*t[7]*t[8]   - t[3]*t[4]*t[9]
                  + t[0]*t[7]*t[9]   + t[1]*t[4]*t[11]  - t[0]*t[5]*t[11];
    cofactor[12] =  t[6]*t[9]*t[12]  - t[5]*t[10]*t[12] - t[6]*t[8]*t[13]
                  + t[4]*t[10]*t[13] + t[5]*t[8]*t[14]  - t[4]*t[9]*t[14];
    cofactor[13] =  t[1]*t[10]*t[12] - t[2]*t[9]*t[12]  + t[2]*t[8]*t[13]
                  - t[0]*t[10]*t[13] - t[1]*t[8]*t[14]  + t[0]*t[9]*t[14];
    cofactor[14] =  t[2]*t[5]*t[12]  - t[1]*t[6]*t[12]  - t[2]*t[4]*t[13]
                  + t[0]*t[6]*t[13]  + t[1]*t[4]*t[14]  - t[0]*t[5]*t[14];
    cofactor[15] =  t[1]*t[6]*t[8]   - t[2]*t[5]*t[8]   + t[2]*t[4]*t[9]
                  - t[0]*t[6]*t[9]   - t[1]*t[4]*t[10]  + t[0]*t[5]*t[10];

    return cofactor/det(t);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
