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

\*---------------------------------------------------------------------------*/

#ifndef symmTensorBlockLduMatrix_H
#define symmTensorBlockLduMatrix_H

#include "coeffFields.H"
#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void Foam::BlockLduMatrix<Foam::symmTensor>::sumDiag()
{
    // Decoupled version
    decoupledSumDiag();
}


template<>
void Foam::BlockLduMatrix<Foam::symmTensor>::negSumDiag()
{
    // Decoupled version
    decoupledNegSumDiag();
}


template<>
void Foam::BlockLduMatrix<symmTensor>::check() const
{
    // Decoupled version
    decoupledCheck();
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
    decoupledRelax(x, b, alpha);
}


template<>
void Foam::BlockLduMatrix<Foam::symmTensor>::operator*=(const scalarField& sf)
{
    // Decoupled version
    decoupledMultEqOp(sf);
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
