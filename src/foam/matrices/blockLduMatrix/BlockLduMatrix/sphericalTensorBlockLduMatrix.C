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
    decoupledSumDiag();
}


template<>
void Foam::BlockLduMatrix<Foam::sphericalTensor>::negSumDiag()
{
    // Decoupled version
    decoupledNegSumDiag();
}


template<>
void Foam::BlockLduMatrix<sphericalTensor>::check() const
{
    // Decoupled version
    decoupledCheck();
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
    decoupledRelax(x, b, alpha);
}


template<>
void Foam::BlockLduMatrix<Foam::sphericalTensor>::operator*=
(
    const scalarField& sf
)
{
    // Decoupled version
    decoupledMultEqOp(sf);
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
