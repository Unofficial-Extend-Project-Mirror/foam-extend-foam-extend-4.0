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
    Template specialisation for scalar block matrix

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    scalarBlockLduMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef scalarBlockLduMatrix_H
#define scalarBlockLduMatrix_H

#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void BlockLduMatrix<scalar>::sumDiag();

template<>
void BlockLduMatrix<scalar>::negSumDiag();

template<>
void BlockLduMatrix<scalar>::check() const;

template<>
void BlockLduMatrix<scalar>::relax
(
    const scalarField& x,
    scalarField& b,
    const scalar alpha
);

template<>
void BlockLduMatrix<scalar>::operator*=(const scalarField& sf);

template<>
void BlockLduMatrix<scalar>::AmulCore
(
    scalarField& mul,
    const scalarField& x
) const;

template<>
void BlockLduMatrix<scalar>::TmulCore
(
    scalarField& mul,
    const scalarField& x
) const;

template<>
void BlockLduMatrix<scalar>::segregateB
(
    scalarField& mul,
    const scalarField& x
) const;

template<>
tmp<scalarField> BlockLduMatrix<scalar>::H(const scalarField& x) const;

template<>
tmp<scalarField> BlockLduMatrix<scalar>::faceH(const scalarField& x) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
