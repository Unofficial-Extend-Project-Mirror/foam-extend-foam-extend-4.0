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

#ifndef sphericalTensorBlockLduMatrix_H
#define sphericalTensorBlockLduMatrix_H

#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void Foam::BlockLduMatrix<Foam::sphericalTensor>::sumDiag()
{
    // Decoupled version
    this->decoupledSumDiag();
}


template<>
void Foam::BlockLduMatrix<Foam::sphericalTensor>::negSumDiag()
{
    // Decoupled version
    this->decoupledNegSumDiag();
}


template<>
void Foam::BlockLduMatrix<sphericalTensor>::check() const
{
    // Decoupled version
    this->decoupledCheck();
}


template<>
void Foam::BlockLduMatrix<Foam::sphericalTensor>::relax
(
    const sphericalTensorField& x,
    sphericalTensorField& b,
    const scalar alpha
)
{
    // Decoupled version
    this->decoupledRelax(x, b, alpha);
}


template<>
void Foam::BlockLduMatrix<Foam::sphericalTensor>::operator*=
(
    const scalarField& sf
)
{
    // Decoupled version
    this->decoupledMultEqOp(sf);
}


template<>
void Foam::BlockLduMatrix<Foam::sphericalTensor>::AmulCore
(
    sphericalTensorField& Ax,
    const sphericalTensorField& x
) const
{
    decoupledAmulCore(Ax, x);
}


template<>
void Foam::BlockLduMatrix<Foam::sphericalTensor>::TmulCore
(
    sphericalTensorField& Tx,
    const sphericalTensorField& x
) const
{
    // Decoupled version
    decoupledTmulCore(Tx, x);
}


template<>
void Foam::BlockLduMatrix<Foam::sphericalTensor>::segregateB
(
    sphericalTensorField&,
    const sphericalTensorField&
) const
{
    FatalErrorIn
    (
        "void Foam::BlockLduMatrix<sphericalTensor>::segregateB\n"
        "(\n"
        "    sphericalTensorField&,\n"
        "    const sphericalTensorField&\n"
        ") const"
    )   << "Requested decoupling of sphericalTensor matrix - never coupled"
        << abort(FatalError);
}


template<>
Foam::tmp<Foam::sphericalTensorField>
Foam::BlockLduMatrix<Foam::sphericalTensor>::H
(
    const sphericalTensorField& x
) const
{
    // Decoupled version
    return decoupledH(x);
}


template<>
Foam::tmp<Foam::sphericalTensorField>
Foam::BlockLduMatrix<Foam::sphericalTensor>::faceH
(
    const sphericalTensorField& x
) const
{
    // Decoupled version
    return decoupledFaceH(x);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
