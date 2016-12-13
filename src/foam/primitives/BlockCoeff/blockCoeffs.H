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

Class
    blockCoeffs

Description

\*---------------------------------------------------------------------------*/

#ifndef blockCoeffs_H
#define blockCoeffs_H

#include "scalarBlockCoeff.H"
#include "sphericalTensorBlockCoeff.H"
#include "symmTensorBlockCoeff.H"
#include "tensorBlockCoeff.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef BlockCoeff<scalar> scalarBlockCoeff;
typedef BlockCoeff<vector> vectorBlockCoeff;
typedef BlockCoeff<sphericalTensor> sphericalTensorBlockCoeff;
typedef BlockCoeff<symmTensor> symmTensorBlockCoeff;
typedef BlockCoeff<tensor> tensorBlockCoeff;


/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

template<class Type>
scalar sumMult
(
    const Field<Type>& f1,
    const Field<Type>& f2
)
{
    Type result = pTraits<Type>::zero;

    typename BlockCoeff<Type>::multiply mult;

    for (register label i = 0; i < f1.size(); i++)
    {
        result += mult(f1[i], f2[i]);
    }

    return cmptSum(result);
}


template<class Type>
scalar gSumMult
(
    const Field<Type>& f1,
    const Field<Type>& f2
)
{
    scalar SumProd = sumMult(f1, f2);
    reduce(SumProd, sumOp<scalar>());
    return SumProd;
}


template<>
inline scalar sumMult
(
    const scalarField& f1,
    const scalarField& f2
)
{
    return sumProd(f1, f2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
