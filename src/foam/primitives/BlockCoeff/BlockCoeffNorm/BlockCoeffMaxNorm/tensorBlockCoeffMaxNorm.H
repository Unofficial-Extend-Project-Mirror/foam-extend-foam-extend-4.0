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
    BlockCoeffMaxNorm

Description
    Implemented to avoid issues with asScalar, asSquare etc.
    This is basically needed since there are specializations of the BlockCoeff
    class.

Author
    Klas Jareteg, 2013-03-07

SourceFiles
    BlockCoeffMaxNorm.C

\*---------------------------------------------------------------------------*/

#ifndef tensorBlockCoeffMaxNorm_H
#define tensorBlockCoeffMaxNorm_H

#include "blockCoeffs.H"
#include "blockCoeffNorms.H"
#include "BlockCoeffNorm.H"
#include "BlockCoeffMaxNorm.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class BlockCoeffMaxNorm Declaration
\*---------------------------------------------------------------------------*/

template<>
inline scalar BlockCoeffMaxNorm<tensor>::normalize
(
    const BlockCoeff<tensor>& a
)
{
    if (a.activeType() == BlockCoeff<tensor>::SCALAR)
    {
        return mag(a.asScalar());
    }
    else if (a.activeType() == BlockCoeff<tensor>::LINEAR)
    {
        return cmptMax(cmptMag(a.asLinear()));
    }
    else
    {
        FatalErrorIn
        (
            "scalar BlockCoeffMaxNorm<tensor>(const BlockCoeff<tensor>& a)"
        )   << "Unknown type" << abort(FatalError);

        return 0;
    }
}


template<>
inline void BlockCoeffMaxNorm<tensor>::coeffMag
(
    const CoeffField<tensor>& a,
    Field<scalar>& b
)
{
    if (a.activeType() == BlockCoeff<tensor>::SCALAR)
    {
        b = mag(a.asScalar());
    }
    else if (a.activeType() == BlockCoeff<tensor>::LINEAR)
    {
        b = cmptMax(cmptMag(a.asLinear()));
    }
    else
    {
        FatalErrorIn
        (
            "scalar BlockCoeffMaxNorm<Type>(const BlockCoeff<tensor>& b)"
        )   << "Unknown type" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
