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
    Specilization for scalar. Implemented to avoid issues with asScalar,
    asSquare etc.
    This is basically needed since there are specializations of the BlockCoeff
    class.

Author
    Klas Jareteg, 2013-03-07

SourceFiles
    BlockCoeffMaxNorm.C

\*---------------------------------------------------------------------------*/

#ifndef scalarBlockCoeffMaxNorm_H
#define scalarBlockCoeffMaxNorm_H

#include "BlockCoeff.H"
#include "CoeffField.H"
#include "BlockCoeffNorm.H"
#include "BlockCoeffMaxNorm.H"
#include "runTimeSelectionTables.H"
#include "scalarCoeffField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class BlockCoeffMaxNorm Declaration
\*---------------------------------------------------------------------------*/

template<>
inline scalar BlockCoeffMaxNorm<scalar>::normalize
(
    const BlockCoeff<scalar>& a
)
{
    return mag(a.asScalar());
}


template<>
inline void BlockCoeffMaxNorm<scalar>::coeffMag
(
    const CoeffField<scalar>& a,
    Field<scalar>& b
)
{
    b = mag(a.asScalar());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
