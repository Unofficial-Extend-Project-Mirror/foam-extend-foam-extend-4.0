/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    BlockDiagGaussSeidelPrecon

Description
    Template specialisation for scalar block diagonally-corrected
    Gauss-Seidel preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef scalarBlockDiagGaussSeidelPrecon_H
#define scalarBlockDiagGaussSeidelPrecon_H

#include "BlockDiagGaussSeidelPrecon.H"
#include "scalarBlockDiagGaussSeidelPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void BlockDiagGaussSeidelPrecon<scalar>::calcInvDiag()
{
    const unallocLabelList& l = this->matrix_.lduAddr().lowerAddr();
    const unallocLabelList& u = this->matrix_.lduAddr().upperAddr();

    // Get diagonal
    const scalarField& d = matrix_.diag();

    scalarField sumMagOffDiag(d.size(), 0);

    // Sum up off-diagonal part of the matrix

    if (matrix_.symmetric())
    {
        // Symmetric matrix

        const scalarField& upper = matrix_.upper();

        for (register label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            sumMagOffDiag[l[coeffI]] += upper[coeffI];
            sumMagOffDiag[u[coeffI]] += upper[coeffI];
        }
    }
    else if (matrix_.asymmetric())
    {
        const scalarField& lower = matrix_.lower();
        const scalarField& upper = matrix_.upper();
    
        for (register label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            sumMagOffDiag[l[coeffI]] += lower[coeffI];
            sumMagOffDiag[u[coeffI]] += upper[coeffI];
        }
    }

    // Compare sum off-diag with diagonal and take larger for the inverse
    // The difference between the two is stored in LUDiag_ for correction

    scalarField signActiveDiag = sign(d);
    scalarField magActiveDiag = mag(d);

    scalarField magNewDiag = Foam::max(magActiveDiag, sumMagOffDiag);

    // Direct inversion of diagonal is sufficient, as the diagonal
    // is linear.  HJ, 20/Aug/2015
    invDiag_ = signActiveDiag/magActiveDiag;

    // Store correction into LUDiag_
    LUDiag_ = signActiveDiag*Foam::max(magNewDiag - magActiveDiag, scalar(0));
}


template<>
void BlockDiagGaussSeidelPrecon<scalar>::precondition
(
    scalarField& x,
    const scalarField& b
) const
{
    // Prepare correction
    if (bPlusLU_.empty())
    {
        bPlusLU_.setSize(b.size(), 0);
    }

    // Add diag coupling to b

    // Multiply overwrites bPlusLU_: no need to initialise
    // Change of sign accounted via change of sign of bPlusLU_
    // HJ, 10/Jan/2016
    multiply(bPlusLU_, LUDiag_, x);
    bPlusLU_ += b;

    if (matrix_.diagonal())
    {
        x = bPlusLU_*invDiag_;
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
            bPlusLU_
        );
    }
}


template<>
void BlockDiagGaussSeidelPrecon<scalar>::preconditionT
(
    scalarField& xT,
    const scalarField& bT
) const
{
    // Prepare correction
    if (bPlusLU_.empty())
    {
        bPlusLU_.setSize(bT.size(), 0);
    }

    // Add diag coupling to b

    // Multiply overwrites bPlusLU_: no need to initialise
    // Change of sign accounted via change of sign of bPlusLU_
    // HJ, 10/Jan/2016
    multiply(bPlusLU_, LUDiag_, xT);
    bPlusLU_ += bT;

    if (matrix_.diagonal())
    {
        xT = bPlusLU_*invDiag_;
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
            bPlusLU_
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
