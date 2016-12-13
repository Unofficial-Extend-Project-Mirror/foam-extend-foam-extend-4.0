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
    ILUC0

Description
    ILU preconditioning without fill in based on Crout algorithm. L and U are
    calculated and stored.

    Reference: Saad, Y.: Iterative Methods for Sparse Linear Systems (2nd
    Edition), SIAM, 2003.

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved

\*---------------------------------------------------------------------------*/

#include "ILUC0.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ILUC0, 0);

    lduPreconditioner::
        addasymMatrixConstructorToTable<ILUC0>
        addILUC0ditionerAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::ILUC0::calcFactorization()
{
    if (!matrix_.diagonal())
    {
        // Get necessary const access to matrix addressing
        const lduAddressing& addr = matrix_.lduAddr();

        // Get upper/lower addressing
        const label* const __restrict__ uPtr = addr.upperAddr().begin();
        const label* const __restrict__ lPtr = addr.lowerAddr().begin();

        // Get owner start addressing
        const label* const __restrict__ ownStartPtr =
            addr.ownerStartAddr().begin();

        // Get losort and losort start addressing
        const label* const __restrict__ lsrPtr = addr.losortAddr().begin();
        const label* const __restrict__ lsrStartPtr =
            addr.losortStartAddr().begin();

        // Get access to factored matrix entries
        scalar* __restrict__ diagPtr = preconDiag_.begin();
        scalar* __restrict__ upperPtr = preconUpper_.begin();
        scalar* __restrict__ lowerPtr = preconLower_.begin();

        // Get access to working fields
        scalar* __restrict__ zPtr = z_.begin();
        scalar* __restrict__ wPtr = w_.begin();

        // Get number of rows
        const label nRows = preconDiag_.size();

        // Define start and end face of this row/column, and number of non zero
        // off diagonal entries
        register label fStart, fEnd, fLsrStart, fLsrEnd;

        // Crout LU factorization

        // Row by row loop (k - loop).
        for (register label rowI = 0; rowI < nRows; ++rowI)
        {
            // Start and end of k-th row (upper) and k-th column (lower)
            fStart = ownStartPtr[rowI];
            fEnd = ownStartPtr[rowI + 1];

            // Initialize temporary working diagonal
            zDiag_ = diagPtr[rowI];

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

            // Lower coeff loop (first i - loop)
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
                zDiag_ -= lowerPtr[losortCoeff]*upperPtr[losortCoeff];

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
            diagRowI = 1.0/zDiag_;

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
            zDiag_ = 0;

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
        forAll (preconDiag_, i)
        {
            preconDiag_[i] = 1.0/preconDiag_[i];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ILUC0::ILUC0
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& dict
)
:
    lduPreconditioner
    (
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces
    ),
    preconDiag_(matrix_.diag()),
    preconLower_(matrix.lower()),
    preconUpper_(matrix.upper()),
    zDiag_(0),
    z_(preconDiag_.size(), 0),
    w_(preconDiag_.size(), 0)
{
    calcFactorization();
}


Foam::ILUC0::ILUC0
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    lduPreconditioner
    (
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces
    ),
    preconDiag_(matrix_.diag()),
    preconLower_(matrix.lower()),
    preconUpper_(matrix.upper()),
    zDiag_(0),
    z_(preconDiag_.size(), 0),
    w_(preconDiag_.size(), 0)
{
    calcFactorization();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ILUC0::~ILUC0()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ILUC0::precondition
(
    scalarField& x,
    const scalarField& b,
    const direction
) const
{
    if (!matrix_.diagonal())
    {
        // Get matrix addressing
        const lduAddressing& addr = matrix_.lduAddr();
        const unallocLabelList& upperAddr = addr.upperAddr();
        const unallocLabelList& lowerAddr = addr.lowerAddr();
        const unallocLabelList& losortAddr = addr.losortAddr();

        // Solve Lz = b with forward substitution. preconLower_ is chosen to
        // be unit triangular. z does not need to be stored

        // Initialize x field
        x = b;

        register label losortCoeffI;
        register label rowI;

        // Forward substitution loop
        forAll (preconLower_, coeffI)
        {
            // Get current losortCoeff to ensure row by row access
            losortCoeffI = losortAddr[coeffI];

            // Subtract already updated lower part from the solution
            x[upperAddr[losortCoeffI]] -=
                preconLower_[losortCoeffI]*x[lowerAddr[losortCoeffI]];
        }

        // Solve Ux = b with back substitution. U is chosen to be upper
        // triangular with diagonal entries corresponding to preconDiag_

        // Multiply with inverse diagonal
        x *= preconDiag_;

        // Back substitution loop
        forAllReverse (preconUpper_, coeffI)
        {
            // Get row index
            rowI = lowerAddr[coeffI];

            // Subtract already updated upper part from the solution
            x[rowI] -=
                preconUpper_[coeffI]*x[upperAddr[coeffI]]*preconDiag_[rowI];
        }
    }
    else
    {
        WarningIn
        (
            "void ILUC0::precondition"
            "(scalarField& x, const scalarField& b, const direction cmpt)"
        )   << "Unnecessary use of ILUC0 preconditioner for diagonal matrix. "
            << nl
            << "Use diagonal preconditioner instead."
            << endl;

        // Diagonal preconditioning
        forAll(x, i)
        {
            x[i] = b[i]*preconDiag_[i];
        }
    }
}


void Foam::ILUC0::preconditionT
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    if (!matrix_.diagonal())
    {
        // Get matrix addressing
        const lduAddressing& addr = matrix_.lduAddr();
        const unallocLabelList& upperAddr = addr.upperAddr();
        const unallocLabelList& lowerAddr = addr.lowerAddr();
        const unallocLabelList& losortAddr = addr.losortAddr();

        // Solve U^T z = b with forward substitution. preconLower_ is chosen to
        // be unit triangular - U^T (transpose U) "contains" diagonal entries. z
        // does not need to be stored.

        // Initialize x field
        forAll(x, i)
        {
            x[i] = b[i]*preconDiag_[i];
        }

        register label losortCoeffI;
        register label rowI;

        // Forward substitution loop
        forAll (preconUpper_, coeffI)
        {
            // Get current losortCoeff to ensure row by row access
            losortCoeffI = losortAddr[coeffI];

            // Get row index
            rowI = upperAddr[losortCoeffI];

            // Subtract already updated lower (upper transpose) part from the
            // solution
            x[rowI] -= preconUpper_[losortCoeffI]*x[lowerAddr[losortCoeffI]]*
                preconDiag_[rowI];
        }

        // Solve L^T x = z with back substitution. L^T is unit upper triangular

        // Back substitution loop
        forAllReverse (preconLower_, coeffI)
        {
            // Subtract already updated upper part from the solution
            x[lowerAddr[coeffI]] -= preconLower_[coeffI]*x[upperAddr[coeffI]];
        }
    }
    else
    {
        WarningIn
        (
            "void ILUC0::preconditionT"
            "(scalarField& x, const scalarField& b, const direction cmpt)"
        )   << "Unnecessary use of ILUC0 preconditioner for diagonal matrix. "
            << nl
            << "Use diagonal preconditioner instead."
            << endl;

        // Diagonal preconditioning
        forAll(x, i)
        {
            x[i] = b[i]*preconDiag_[i];
        }
    }
}


// ************************************************************************* //
