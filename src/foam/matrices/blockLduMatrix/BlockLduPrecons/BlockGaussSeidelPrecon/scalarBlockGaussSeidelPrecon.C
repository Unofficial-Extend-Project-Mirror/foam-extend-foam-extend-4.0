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
    BlockGaussSeidelPrecon

Description
    Template specialisation for scalar block Gauss-Seidel preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef scalarBlockGaussSeidelPrecon_H
#define scalarBlockGaussSeidelPrecon_H

#include "BlockGaussSeidelPrecon.H"
#include "scalarBlockGaussSeidelPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void BlockGaussSeidelPrecon<scalar>::calcInvDiag()
{
    // Direct inversion of diagonal is sufficient, as the diagonal
    // is linear.  HJ, 20/Aug/2015
    invDiag_ = 1/this->matrix_.diag();
}


template<>
void BlockGaussSeidelPrecon<scalar>::precondition
(
    scalarField& x,
    const scalarField& b
) const
{
    if (matrix_.diagonal())
    {
        x = b*invDiag_;
    }
    else if (matrix_.symmetric() || matrix_.asymmetric())
    {
        const scalarField& LowerCoeff = matrix_.lower();
        const scalarField& UpperCoeff = matrix_.upper();

        BlockSweep
        (
            x,
            invDiag_,
            LowerCoeff,
            UpperCoeff,
            b
        );
    }
}


template<>
void BlockGaussSeidelPrecon<scalar>::preconditionT
(
    scalarField& xT,
    const scalarField& bT
) const
{
    if (matrix_.diagonal())
    {
        xT = bT*invDiag_;
    }
    else if (matrix_.symmetric() || matrix_.asymmetric())
    {
        const scalarField& LowerCoeff = matrix_.lower();
        const scalarField& UpperCoeff = matrix_.upper();

        // Swap lower and upper coefficients, transposed matrix
        BlockSweep
        (
            xT,
            invDiag_,
            UpperCoeff,
            LowerCoeff,
            bT
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
