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
    Template specialisation for sphericalTensor block matrix

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    sphericalTensorBlockLduMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalTensorBlockLduMatrix_H
#define sphericalTensorBlockLduMatrix_H

#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void BlockLduMatrix<sphericalTensor>::sumDiag();

template<>
void BlockLduMatrix<sphericalTensor>::negSumDiag();

template<>
void BlockLduMatrix<sphericalTensor>::check() const;

template<>
void BlockLduMatrix<sphericalTensor>::relax
(
    const sphericalTensorField& x,
    sphericalTensorField& b,
    const scalar alpha
);

template<>
void BlockLduMatrix<sphericalTensor>::operator*=(const scalarField& sf);

template<>
void BlockLduMatrix<sphericalTensor>::AmulCore
(
    sphericalTensorField& mul,
    const sphericalTensorField& x
) const;

template<>
void BlockLduMatrix<sphericalTensor>::TmulCore
(
    sphericalTensorField& mul,
    const sphericalTensorField& x
) const;

template<>
void BlockLduMatrix<sphericalTensor>::segregateB
(
    sphericalTensorField& mul,
    const sphericalTensorField& x
) const;


template<>
tmp<sphericalTensorField>
BlockLduMatrix<sphericalTensor>::H(const sphericalTensorField& x) const;

template<>
tmp<sphericalTensorField>
BlockLduMatrix<sphericalTensor>::faceH(const sphericalTensorField& x) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
