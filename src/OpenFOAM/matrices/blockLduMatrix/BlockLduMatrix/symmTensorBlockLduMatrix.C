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

\*---------------------------------------------------------------------------*/

#ifndef symmTensorBlockLduMatrix_H
#define symmTensorBlockLduMatrix_H

#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void Foam::BlockLduMatrix<Foam::symmTensor>::sumDiag()
{
    // Decoupled version
    this->decoupledSumDiag();
}


template<>
void Foam::BlockLduMatrix<Foam::symmTensor>::negSumDiag()
{
    // Decoupled version
    this->decoupledNegSumDiag();
}


template<>
void Foam::BlockLduMatrix<symmTensor>::check() const
{
    // Decoupled version
    this->decoupledCheck();
}


template<>
void Foam::BlockLduMatrix<Foam::symmTensor>::relax
(
    const symmTensorField& x,
    symmTensorField& b,
    const scalar alpha
)
{
    // Decoupled version
    this->decoupledRelax(x, b, alpha);
}


template<>
void Foam::BlockLduMatrix<Foam::symmTensor>::operator*=(const scalarField& sf)
{
    // Decoupled version
    this->decoupledMultEqOp(sf);
}


template<>
void Foam::BlockLduMatrix<Foam::symmTensor>::AmulCore
(
    symmTensorField& Ax,
    const symmTensorField& x
) const
{
    decoupledAmulCore(Ax, x);
}


template<>
void Foam::BlockLduMatrix<Foam::symmTensor>::TmulCore
(
    symmTensorField& Tx,
    const symmTensorField& x
) const
{
    // Decoupled version
    decoupledTmulCore(Tx, x);
}


template<>
void Foam::BlockLduMatrix<Foam::symmTensor>::segregateB
(
    symmTensorField&,
    const symmTensorField&
) const
{
    FatalErrorIn
    (
        "void Foam::BlockLduMatrix<symmTensor>::segregateB\n"
        "(\n"
        "    symmTensorField&,\n"
        "    const symmTensorField&\n"
        ") const"
    )   << "Requested decoupling of symmTensor matrix - never coupled"
        << abort(FatalError);
}


template<>
Foam::tmp<Foam::symmTensorField>
Foam::BlockLduMatrix<Foam::symmTensor>::H(const symmTensorField& x) const
{
    // Decoupled version
    return decoupledH(x);
}


template<>
Foam::tmp<Foam::symmTensorField>
Foam::BlockLduMatrix<Foam::symmTensor>::faceH(const symmTensorField& x) const
{
    // Decoupled version
    return decoupledFaceH(x);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
