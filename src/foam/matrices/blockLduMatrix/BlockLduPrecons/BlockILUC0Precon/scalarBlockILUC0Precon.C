/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    BlockILUC0Precon

Description
    Template specialisation for scalar block ILUCp preconditioning

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef scalarBlockILUC0Precon_H
#define scalarBlockILUC0Precon_H

#include "BlockILUC0Precon.H"
#include "scalarBlockILUC0Precon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
template<>
void BlockILUC0Precon<scalar>::calcActiveTypeFactorization
(
    scalarField& preconD,
    scalarField& preconUpper,
    scalarField& preconLower
) const
{
    if (!matrix_.diagonal())
    {
        // Get number of rows
        const label nRows = preconD.size();

        // Allocate working fields
        scalarField z(nRows, 0.0);
        scalarField w(nRows, 0.0);
        scalar zDiag = 0.0;

        // Get necessary const access to extended ldu addressing
        const lduAddressing& addr = this->matrix_.lduAddr();

        // Get upper/lower extended addressing
        const label* const __restrict__ uPtr = addr.upperAddr().begin();
        const label* const __restrict__ lPtr = addr.lowerAddr().begin();

        // Get extended owner start addressing
        const label* const __restrict__ ownStartPtr =
            addr.ownerStartAddr().begin();

        // Get extended losort and losort start addressing
        const label* const __restrict__ lsrPtr =
            addr.losortAddr().begin();
        const label* const __restrict__ lsrStartPtr =
            addr.losortStartAddr().begin();

        // Get access to factored matrix entries
        scalar* __restrict__ diagPtr = preconD.begin();
        scalar* __restrict__ upperPtr = preconUpper.begin();
        scalar* __restrict__ lowerPtr = preconLower.begin();

        // Get access to working fields
        scalar* __restrict__ zPtr = z.begin();
        scalar* __restrict__ wPtr = w.begin();

        // Define start and end face ("virtual" face when extended addressing is
        // used) of this row/column.
        label fStart, fEnd, fLsrStart, fLsrEnd;

        // Crout LU factorization

        // Row by row loop (k - loop).
        for (label rowI = 0; rowI < nRows; ++rowI)
        {
            // Start and end of k-th row (upper) and k-th column (lower)
            fStart = ownStartPtr[rowI];
            fEnd = ownStartPtr[rowI + 1];

            // Initialize temporary working diagonal
            zDiag = diagPtr[rowI];

            // Initialize temporary working row field
            for (label faceI = fStart; faceI < fEnd; ++faceI)
            {
                // Note: z addressed by neighbour of face (column index for
                // upper), w addressed by neighbour of face (row index for
                // lower)
                zPtr[uPtr[faceI]] = upperPtr[faceI];
                wPtr[uPtr[faceI]] = lowerPtr[faceI];
            }

            // Start and end of k-th row (lower) and k-th column (upper)
            fLsrStart = lsrStartPtr[rowI];
            fLsrEnd = lsrStartPtr[rowI + 1];

            // Lower/upper coeff loop (i - loop)
            for
            (
                label faceLsrI = fLsrStart;
                faceLsrI < fLsrEnd;
                ++faceLsrI
            )
            {
                // Get losort coefficient for this face
                const label losortCoeff = lsrPtr[faceLsrI];

                // Get corresponding row index for upper (i label)
                const label i = lPtr[losortCoeff];

                // Update diagonal
                zDiag -= lowerPtr[losortCoeff]*upperPtr[losortCoeff];

                // Get end of row for cell i
                const label fEndRowi = ownStartPtr[i + 1];

                // Upper coeff loop (additional loop to avoid checking the
                // existence of certain upper coeffs)
                for
                (
                    // Diagonal is already updated (losortCoeff + 1 = start)
                    label faceI = losortCoeff + 1;
                    faceI < fEndRowi;
                    ++faceI
                )
                {
                    zPtr[uPtr[faceI]] -= lowerPtr[losortCoeff]*upperPtr[faceI];
                    wPtr[uPtr[faceI]] -= upperPtr[losortCoeff]*lowerPtr[faceI];
                }
            }

            // Update diagonal entry, inverting it for future use
            scalar& diagRowI = diagPtr[rowI];
            diagRowI = 1.0/zDiag;

            // Index for updating L and U
            label zwIndex;

            // Update upper and lower coeffs
            for (label faceI = fStart; faceI < fEnd; ++faceI)
            {
                // Get index for current face
                zwIndex = uPtr[faceI];

                // Update L and U decomposition for this row (column)
                upperPtr[faceI] = zPtr[zwIndex];
                lowerPtr[faceI] = wPtr[zwIndex]*diagRowI;
            }

            // Reset temporary working fields
            zDiag = 0;

            // Only reset parts of the working fields that have been updated in
            // this step (for this row and column)
            for
            (
                label faceLsrI = fLsrStart;
                faceLsrI < fLsrEnd;
                ++faceLsrI
            )
            {
                // Get losort coefficient for this face
                const label losortCoeff = lsrPtr[faceLsrI];

                // Get corresponding row index for upper (i label)
                const label i = lPtr[losortCoeff];

                // Get end of row for cell i
                const label fEndRowi = ownStartPtr[i + 1];

                for
                (
                    label faceI = losortCoeff + 1;
                    faceI < fEndRowi;
                    ++faceI
                )
                {
                    zPtr[uPtr[faceI]] = 0.0;
                    wPtr[uPtr[faceI]] = 0.0;
                }
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "template<>\n"
            "template<>\n"
            "void BlockILUC0Precon<scalar>::calcFactorization\n"
            "(\n"
            "    scalarField& preconD,\n"
            "    scalarField& extUpper,\n"
            "    scalarField& extLower,\n"
            "    scalarField& zDiag\n,"
            "    scalarField& z,\n"
            "    scalarField& w,\n"
            ") const"
        )   << "Unnecessary use of BlockILUC0 preconditioner for diagonal "
            << "matrix."
            << nl
            << "Use BlockDiagonal preconditioner instead."
            << abort(FatalError);
    }
}


template<>
void BlockILUC0Precon<scalar>::calcFactorization() const
{
    calcActiveTypeFactorization
    (
        preconDiag_.asScalar(),
        preconUpper_.asScalar(),
        preconLower_.asScalar()
    );
}


template<>
void BlockILUC0Precon<scalar>::precondition
(
    scalarField& x,
    const scalarField& b
) const
{
    // Decoupled version
    notImplemented("void Foam::BlockILUC0Precon<scalar>::precondition");
}


template<>
void BlockILUC0Precon<scalar>::preconditionT
(
    scalarField& xT,
    const scalarField& bT
) const
{
    // Decoupled version
    notImplemented("void Foam::BlockILUC0Precon<scalar>::preconditionT");
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
