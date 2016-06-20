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
    Template specialisation for symmTensor block matrix

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    symmTensorBlockLduMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef symmTensorBlockLduMatrix_H
#define symmTensorBlockLduMatrix_H

#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void BlockLduMatrix<symmTensor>::sumDiag();

template<>
void BlockLduMatrix<symmTensor>::negSumDiag();

template<>
void BlockLduMatrix<symmTensor>::check() const;

template<>
void BlockLduMatrix<symmTensor>::relax
(
    const symmTensorField& x,
    symmTensorField& b,
    const scalar alpha
);

template<>
void BlockLduMatrix<symmTensor>::operator*=(const scalarField& sf);

template<>
void BlockLduMatrix<symmTensor>::AmulCore
(
    symmTensorField& mul,
    const symmTensorField& x
) const;

template<>
void BlockLduMatrix<symmTensor>::TmulCore
(
    symmTensorField& mul,
    const symmTensorField& x
) const;

template<>
void BlockLduMatrix<symmTensor>::segregateB
(
    symmTensorField& mul,
    const symmTensorField& x
) const;


template<>
tmp<symmTensorField>
BlockLduMatrix<symmTensor>::H(const symmTensorField& x) const;

template<>
tmp<symmTensorField>
BlockLduMatrix<symmTensor>::faceH(const symmTensorField& x) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
