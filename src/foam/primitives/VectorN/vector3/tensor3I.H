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
    tensor3

Description
    TensorN of 3 scalars.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


//- Return the determinant of a tensor
inline scalar det(const tensor3& t)
{
    return
    (
        t[0]*t[4]*t[8] + t[1]*t[5]*t[6]
      + t[2]*t[3]*t[7] - t[0]*t[5]*t[7]
      - t[1]*t[3]*t[8] - t[2]*t[4]*t[6]
    );
}


//- Return the inverse of a tensor given the determinant
inline tensor3 inv(const tensor3& t)
{
    tensor3 cofactor;

    cofactor[0] = t[4]*t[8] - t[7]*t[5];
    cofactor[1] = t[6]*t[5] - t[3]*t[8];
    cofactor[2] = t[3]*t[7] - t[4]*t[6];
    cofactor[3] = t[2]*t[7] - t[1]*t[8];
    cofactor[4] = t[0]*t[8] - t[2]*t[6];
    cofactor[5] = t[1]*t[6] - t[0]*t[7];
    cofactor[6] = t[1]*t[5] - t[2]*t[4];
    cofactor[7] = t[3]*t[2] - t[0]*t[5];
    cofactor[8] = t[0]*t[4] - t[3]*t[1];

    return cofactor/det(t);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
