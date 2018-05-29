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
    Class for two norm of block coeffs. Specilization for tensor. Implemented
    to avoid issues with asScalar, asSquare etc.
    This is basically needed since there are specializations of the BlockCoeff
    class.

Author
    Klas Jareteg, 2013-01-30

\*---------------------------------------------------------------------------*/

#ifndef tensorBlockCoeffTwoNorm_H
#define tensorBlockCoeffTwoNorm_H

#include "blockCoeffs.H"
#include "blockCoeffNorms.H"
#include "BlockCoeffNorm.H"
#include "BlockCoeffTwoNorm.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
inline scalar BlockCoeffTwoNorm<tensor>::normalize
(
    const BlockCoeff<tensor>& a
)
{
    if (a.activeType() == BlockCoeff<tensor>::SCALAR)
    {
        // Note: for two-norm, use mag
        return sign(a.component(cmpt_))*mag(a.asScalar());
    }
    else if (a.activeType() == BlockCoeff<tensor>::LINEAR)
    {
        // Note: for two-norm, use mag
        return sign(a.component(cmpt_))*mag(a.asLinear());
    }
    else
    {
        FatalErrorIn
        (
            "scalar BlockCoeffTwoNorm<Type>(const BlockCoeff<Type>& a)"
        )   << "Unknown type" << abort(FatalError);

        return 0;
    }
}


template<>
inline void BlockCoeffTwoNorm<tensor>::normalize
(
    Field<scalar>& b,
    const CoeffField<tensor>& a
)
{
    if (a.activeType() == BlockCoeff<tensor>::SCALAR)
    {
        // Note: for two-norm, use mag
        // However, the norm is signed
        // b = mag(a.asScalar());
        b = a.asScalar();
    }
    else if (a.activeType() == BlockCoeff<tensor>::LINEAR)
    {
        // Note: for two-norm, use mag
        b = sign(a.component(cmpt_))*mag(a.asLinear());
    }
    else
    {
        FatalErrorIn
        (
            "scalar BlockCoeffTwoNorm<Type>(const BlockCoeff<Type>& b)"
        )   << "Unknown type" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
