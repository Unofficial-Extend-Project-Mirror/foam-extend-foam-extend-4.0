/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

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
