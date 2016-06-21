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
    BlockLduMatrix

Description
    Template specialisation for tensor block matrix

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    tensorBlockLduMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef tensorBlockLduMatrix_H
#define tensorBlockLduMatrix_H

#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void BlockLduMatrix<tensor>::sumDiag();

template<>
void BlockLduMatrix<tensor>::negSumDiag();

template<>
void BlockLduMatrix<tensor>::check() const;

template<>
void BlockLduMatrix<tensor>::relax
(
    const tensorField& x,
    tensorField& b,
    const scalar alpha
);

template<>
void BlockLduMatrix<tensor>::operator*=(const scalarField& sf);

template<>
void BlockLduMatrix<tensor>::AmulCore
(
    tensorField& mul,
    const tensorField& x
) const;

template<>
void BlockLduMatrix<tensor>::TmulCore
(
    tensorField& mul,
    const tensorField& x
) const;

template<>
void BlockLduMatrix<tensor>::segregateB
(
    tensorField& mul,
    const tensorField& x
) const;


template<>
tmp<tensorField> BlockLduMatrix<tensor>::H(const tensorField& x) const;

template<>
tmp<tensorField> BlockLduMatrix<tensor>::faceH(const tensorField& x) const;



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
