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

#ifndef tensorBlockLduMatrix_H
#define tensorBlockLduMatrix_H

#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void Foam::BlockLduMatrix<Foam::tensor>::sumDiag()
{
    // Decoupled version
    decoupledSumDiag();
}


template<>
void Foam::BlockLduMatrix<Foam::tensor>::negSumDiag()
{
    // Decoupled version
    decoupledNegSumDiag();
}


template<>
void Foam::BlockLduMatrix<tensor>::check() const
{
    // Decoupled version
    decoupledCheck();
}


template<>
void Foam::BlockLduMatrix<Foam::tensor>::relax
(
    const tensorField& x,
    tensorField& b,
    const scalar alpha
)
{
    // Decoupled version
    decoupledRelax(x, b, alpha);
}


template<>
void Foam::BlockLduMatrix<Foam::tensor>::operator*=(const scalarField& sf)
{
    // Decoupled version
    decoupledMultEqOp(sf);
}


template<>
void Foam::BlockLduMatrix<Foam::tensor>::AmulCore
(
    tensorField& Ax,
    const tensorField& x
) const
{
    decoupledAmulCore(Ax, x);
}


template<>
void Foam::BlockLduMatrix<Foam::tensor>::TmulCore
(
    tensorField& Tx,
    const tensorField& x
) const
{
    // Decoupled version
    decoupledTmulCore(Tx, x);
}


template<>
void Foam::BlockLduMatrix<Foam::tensor>::segregateB
(
    tensorField&,
    const tensorField&
) const
{
    FatalErrorIn
    (
        "void Foam::BlockLduMatrix<tensor>::segregateB\n"
        "(\n"
        "    tensorField&,\n"
        "    const tensorField&\n"
        ") const"
    )   << "Requested decoupling of tensor matrix - never coupled"
        << abort(FatalError);
}


template<>
Foam::tmp<Foam::tensorField>
Foam::BlockLduMatrix<Foam::tensor>::H(const tensorField& x) const
{
    // Decoupled version
    return decoupledH(x);
}


template<>
Foam::tmp<Foam::tensorField>
Foam::BlockLduMatrix<Foam::tensor>::faceH(const tensorField& x) const
{
    // Decoupled version
    return decoupledFaceH(x);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
