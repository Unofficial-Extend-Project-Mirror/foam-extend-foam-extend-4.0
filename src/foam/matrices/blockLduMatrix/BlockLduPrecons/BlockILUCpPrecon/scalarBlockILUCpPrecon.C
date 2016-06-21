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
    BlockILUCpPrecon

Description
    Template specialisation for scalar block ILUCp preconditioning

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef scalarBlockILUCpPrecon_H
#define scalarBlockILUCpPrecon_H

#include "BlockILUCpPrecon.H"
#include "scalarBlockILUCpPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
template<>
void BlockILUCpPrecon<scalar>::calcActiveTypeFactorization
(
    scalarField& preconD,
    scalarField& extUpper,
    scalarField& extLower
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
        const extendedLduAddressing& addr = extBlockMatrix_.extendedLduAddr();

        // Get upper/lower extended addressing
        const label* const __restrict__ uPtr = addr.extendedUpperAddr().begin();
        const label* const __restrict__ lPtr = addr.extendedLowerAddr().begin();

        // Get extended owner start addressing
        const label* const __restrict__ ownStartPtr =
            addr.extendedOwnerStartAddr().begin();

        // Get extended losort and losort start addressing
        const label* const __restrict__ lsrPtr =
            addr.extendedLosortAddr().begin();
        const label* const __restrict__ lsrStartPtr =
            addr.extendedLosortStartAddr().begin();

        // Get access to factored matrix entries
        scalar* __restrict__ diagPtr = preconD.begin();
        scalar* __restrict__ upperPtr = extUpper.begin();
        scalar* __restrict__ lowerPtr = extLower.begin();

        // Get access to working fields
        scalar* __restrict__ zPtr = z.begin();
        scalar* __restrict__ wPtr = w.begin();

        // Define start and end face ("virtual" face when extended addressing is
        // used) of this row/column.
        register label fStart, fEnd, fLsrStart, fLsrEnd;

        // Crout LU factorization

        // Row by row loop (k - loop).
        for (register label rowI = 0; rowI < nRows; ++rowI)
        {
            // Start and end of k-th row (upper) and k-th column (lower)
            fStart = ownStartPtr[rowI];
            fEnd = ownStartPtr[rowI + 1];

            // Initialize temporary working diagonal
            zDiag = diagPtr[rowI];

            // Initialize temporary working row field
            for (register label faceI = fStart; faceI < fEnd; ++faceI)
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
                register label faceLsrI = fLsrStart;
                faceLsrI < fLsrEnd;
                ++faceLsrI
            )
            {
                // Get losort coefficient for this face
                const register label losortCoeff = lsrPtr[faceLsrI];

                // Get corresponding row index for upper (i label)
                const label i = lPtr[losortCoeff];

                // Update diagonal
                zDiag -= lowerPtr[losortCoeff]*upperPtr[losortCoeff];

                // Get end of row for cell i
                const register label fEndRowi = ownStartPtr[i + 1];

                // Upper coeff loop (additional loop to avoid checking the
                // existence of certain upper coeffs)
                for
                (
                    // Diagonal is already updated (losortCoeff + 1 = start)
                    register label faceI = losortCoeff + 1;
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
            register label zwIndex;

            // Update upper and lower coeffs
            for (register label faceI = fStart; faceI < fEnd; ++faceI)
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
                register label faceLsrI = fLsrStart;
                faceLsrI < fLsrEnd;
                ++faceLsrI
            )
            {
                // Get losort coefficient for this face
                const register label losortCoeff = lsrPtr[faceLsrI];

                // Get corresponding row index for upper (i label)
                const label i = lPtr[losortCoeff];

                // Get end of row for cell i
                const register label fEndRowi = ownStartPtr[i + 1];

                for
                (
                    register label faceI = losortCoeff + 1;
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
            "void BlockILUCpPrecon<scalar>::calcFactorization\n"
            "(\n"
            "    scalarField& preconD,\n"
            "    scalarField& extUpper,\n"
            "    scalarField& extLower,\n"
            "    scalarField& zDiag\n,"
            "    scalarField& z,\n"
            "    scalarField& w,\n"
            ") const"
        )   << "Unnecessary use of BlockILUCp preconditioner for diagonal "
            << "matrix."
            << nl
            << "Use BlockDiagonal preconditioner instead."
            << abort(FatalError);
    }
}


template<>
void BlockILUCpPrecon<scalar>::calcFactorization() const
{
    calcActiveTypeFactorization
    (
        preconDiag_.asScalar(),
        extBlockMatrix_.extendedUpper().asScalar(),
        extBlockMatrix_.extendedLower().asScalar()
    );
}


template<>
void BlockILUCpPrecon<scalar>::precondition
(
    scalarField& x,
    const scalarField& b
) const
{
    // Decoupled version
    notImplemented("void Foam::BlockILUCpPrecon<scalar>::precondition");
}


template<>
void BlockILUCpPrecon<scalar>::preconditionT
(
    scalarField& xT,
    const scalarField& bT
) const
{
    // Decoupled version
    notImplemented("void Foam::BlockILUCpPrecon<scalar>::preconditionT");
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
