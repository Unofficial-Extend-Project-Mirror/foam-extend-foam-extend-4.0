/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    BlockCoeffTwoNorm

Description
    Class for two norm of block coeffs. Specilization for scalar. Implemented
    to avoid issues with asScalar, asSquare etc.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef scalarBlockCoeffTwoNorm_H
#define scalarBlockCoeffTwoNorm_H

#include "blockCoeffs.H"
#include "blockCoeffNorms.H"
#include "BlockCoeffNorm.H"
#include "BlockCoeffTwoNorm.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
inline scalar BlockCoeffTwoNorm<scalar>::normalize
(
    const BlockCoeff<scalar>& a
)
{
    // Note: for two-norm, mag should be used.
    // However, the norm is signed
    // return mag(a.asScalar());
    return a.asScalar();
}


template<>
inline void BlockCoeffTwoNorm<scalar>::normalize
(
    Field<scalar>& b,
    const CoeffField<scalar>& a
)
{
    // Note: for two-norm, mag should be used.
    // However, the norm is signed
    // b = mag(a.asScalar());
    b = a.asScalar();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
