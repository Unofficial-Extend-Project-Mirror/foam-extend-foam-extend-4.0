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

Description
    Block variant of ILU preconditioning with arbitrary level of fill in (p),
    based on Crout algorithm.

    Reference: Saad, Y.: Iterative Methods for Sparse Linear Systems
    (2nd Edition), SIAM, 2003.

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "BlockILUC0Precon.H"
#include "error.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class LDUType>
void Foam::BlockILUC0Precon<Type>::calcActiveTypeFactorization
(
    Field<LDUType>& preconD,
    Field<LDUType>& preconUpper,
    Field<LDUType>& preconLower
) const
{
    if (!this->matrix_.diagonal())
    {
        // Get number of rows
        const label nRows = preconD.size();

        // Allocate working fields
        Field<LDUType> z(nRows, pTraits<LDUType>::zero);
        Field<LDUType> w(nRows, pTraits<LDUType>::zero);
        LDUType zDiag(pTraits<LDUType>::zero);

        // Create multiplication function object
        typename BlockCoeff<Type>::multiply mult;

        // Get necessary const access to matrix addressing
        const lduAddressing& addr = this->matrix_.lduAddr();

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
        LDUType* __restrict__ diagPtr = preconD.begin();
        LDUType* __restrict__ upperPtr = preconUpper.begin();
        LDUType* __restrict__ lowerPtr = preconLower.begin();

        // Get access to working fields
        LDUType* __restrict__ zPtr = z.begin();
        LDUType* __restrict__ wPtr = w.begin();

        // Define start and end face ("virtual" face when extended addressing
        // is used) of this row/column.
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

            // Initialize temporary working row and column fields
            for (label faceI = fStart; faceI < fEnd; ++faceI)
            {
                // Note: z addressed by neighbour of face (column index for
                // upper)
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
                zDiag -= mult.activeTypeMultiply
                (
                    lowerPtr[losortCoeff],
                    upperPtr[losortCoeff]
                );

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
                    zPtr[uPtr[faceI]] -= mult.activeTypeMultiply
                    (
                        lowerPtr[losortCoeff],
                        upperPtr[faceI]
                    );

                    wPtr[uPtr[faceI]] -= mult.activeTypeMultiply
                    (
                        lowerPtr[faceI],
                        upperPtr[losortCoeff]
                    );
                }
            }

            // Update diagonal entry, inverting it for future use
            LDUType& diagRowI = diagPtr[rowI];
            diagRowI = mult.inverse(zDiag);

            // Index for updating L and U
            label zwIndex;

            // Update upper and lower coeffs
            for (label faceI = fStart; faceI < fEnd; ++faceI)
            {
                // Get index for current face
                zwIndex = uPtr[faceI];

                // Update L and U decomposition for this row (column)
                upperPtr[faceI] = zPtr[zwIndex];
                lowerPtr[faceI] = mult.activeTypeMultiply
                (
                    wPtr[zwIndex],
                    diagRowI
                );
            }

            // Reset temporary working fields
            zDiag = pTraits<LDUType>::zero;

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
                    zPtr[uPtr[faceI]] = pTraits<LDUType>::zero;
                    wPtr[uPtr[faceI]] = pTraits<LDUType>::zero;
                }
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "template <class Type>\n"
            "template <class LDUType>\n"
            "void BlockILUC0Precon<Type>::calcActiveTypeFactorization\n"
            "(\n"
            "    Field<LDUType>& preconD,\n"
            "    Field<LDUType>& preconUpper,\n"
            "    Field<LDUType>& preconLower,\n"
            ") const"
        )   << "Unnecessary use of BlockILUC0 preconditioner for diagonal "
            << "matrix."
            << nl
            << "Use BlockDiagonal preconditioner instead."
            << abort(FatalError);
    }
}


template<class Type>
void Foam::BlockILUC0Precon<Type>::calcFactorization() const
{
    // Get upper and lower matrix factors
    CoeffField<Type>& Lower = preconLower_;
    CoeffField<Type>& Upper = preconUpper_;

    // Calculate factorization
    // Note: lower, diag and upper must have same type as required by the
    // algorithm. This is handled by lowest possible promotion
    if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
    {
        if (Upper.activeType() == blockCoeffBase::SCALAR)
        {
            calcActiveTypeFactorization
            (
                preconDiag_.asScalar(),
                Upper.asScalar(),
                Lower.asScalar()
            );
        }
        else if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
            calcActiveTypeFactorization
            (
                preconDiag_.asLinear(), // Promotes to linear
                Upper.asLinear(),
                Lower.asLinear()
            );
        }
        else if (Upper.activeType() == blockCoeffBase::SQUARE)
        {
            calcActiveTypeFactorization
            (
                preconDiag_.asSquare(), // Promotes to square
                Upper.asSquare(),
                Lower.asSquare()
            );
        }
    }
    else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
    {
        if (Upper.activeType() == blockCoeffBase::SCALAR)
        {
            calcActiveTypeFactorization
            (
                preconDiag_.asLinear(),
                Upper.asLinear(), // Promotes to linear
                Lower.asLinear()  // Promotes to linear
            );
        }
        else if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
            calcActiveTypeFactorization
            (
                preconDiag_.asLinear(),
                Upper.asLinear(),
                Lower.asLinear()
            );
        }
        else if (Upper.activeType() == blockCoeffBase::SQUARE)
        {
            calcActiveTypeFactorization
            (
                preconDiag_.asSquare(), // Promotes to square
                Upper.asSquare(),
                Lower.asSquare()
            );
        }
    }
    else if (preconDiag_.activeType() == blockCoeffBase::SQUARE)
    {
        if (Upper.activeType() == blockCoeffBase::SCALAR)
        {
            calcActiveTypeFactorization
            (
                preconDiag_.asSquare(),
                Upper.asSquare(), // Promotes to square
                Lower.asSquare()  // Promotes to square
            );
        }
        else if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
            calcActiveTypeFactorization
            (
                preconDiag_.asSquare(),
                Upper.asSquare(), // Promotes to square
                Lower.asSquare()  // Promotes to square
            );
        }
        else if (Upper.activeType() == blockCoeffBase::SQUARE)
        {
            calcActiveTypeFactorization
            (
                preconDiag_.asSquare(),
                Upper.asSquare(),
                Lower.asSquare()
            );
        }
    }
}


template<class Type>
template<class LDUType>
void Foam::BlockILUC0Precon<Type>::LUSubstitute
(
    Field<Type>& x,
    const Field<LDUType>& preconD,
    const Field<LDUType>& upper,
    const Field<LDUType>& lower,
    const Field<Type>& b
) const
{
    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    // Get matrix addressing
    const lduAddressing& addr = this->matrix_.lduAddr();
    const unallocLabelList& upperAddr = addr.upperAddr();
    const unallocLabelList& lowerAddr = addr.lowerAddr();
    const unallocLabelList& losortAddr = addr.losortAddr();

    // Solve Lz = b with forward substitution in block form. lower is chosen
    // to be unit triangular. z does not need to be stored

    // Initialize x field
    x = b;

    label losortCoeffI;
    label rowI;

    // Forward substitution loop
    forAll (lower, coeffI)
    {
        // Get current losortCoeff to ensure row by row access
        losortCoeffI = losortAddr[coeffI];

        // Subtract already updated lower part from the solution
        x[upperAddr[losortCoeffI]] -= mult
        (
            lower[losortCoeffI],
            x[lowerAddr[losortCoeffI]]
        );
    }

    // Solve Ux = b with back substitution in block form. U is chosen to be
    // upper triangular with diagonal entries corresponding to preconD

    // Multiply with inverse diagonal
    forAll(x, i)
    {
        x[i] = mult(preconD[i], x[i]);
    }

    // Back substitution loop
    forAllReverse (upper, coeffI)
    {
        // Get row index
        rowI = lowerAddr[coeffI];

        // Subtract already updated upper part from the solution
        x[rowI] -= mult
        (
            preconD[rowI],
            mult
            (
                upper[coeffI],
                x[upperAddr[coeffI]]
            )
        );
    }
}


template<class Type>
template<class LDUType>
void Foam::BlockILUC0Precon<Type>::LUSubstituteT
(
    Field<Type>& xT,
    const Field<LDUType>& preconD,
    const Field<LDUType>& upper,
    const Field<LDUType>& lower,
    const Field<Type>& bT
) const
{
    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    // Get matrix addressing
    const lduAddressing& addr = this->matrix_.lduAddr();
    const unallocLabelList& upperAddr = addr.upperAddr();
    const unallocLabelList& lowerAddr = addr.lowerAddr();
    const unallocLabelList& losortAddr = addr.losortAddr();

    // Solve U^T z = b with forward substitution in block form. lower is
    // chosen to be unit triangular - U^T (transpose U) "contains" diagonal
    // entries. z does not need to be stored

    // Note: transpose should be used for all block coeffs.

    // Initialize x field
    forAll(xT, i)
    {
        xT[i] = mult
        (
            mult.transpose(preconD[i]),
            bT[i]
        );
    }

    label losortCoeffI;
    label rowI;

    // Forward substitution loop
    forAll (upper, coeffI)
    {
        // Get current losortCoeff to ensure row by row access
        losortCoeffI = losortAddr[coeffI];

        // Get row index
        rowI = upperAddr[losortCoeffI];

        // Subtract already updated lower (upper transpose) part from the
        // solution
        xT[rowI] -= mult
        (
            mult.transpose(preconD[rowI]),
            mult
            (
                mult.transpose(upper[losortCoeffI]),
                xT[lowerAddr[losortCoeffI]]
            )
        );
    }

    // Solve L^T x = z with back substitution. L^T is unit upper triangular

    // Back substitution loop
    forAllReverse (lower, coeffI)
    {
        // Subtract already updated upper part from the solution
        xT[lowerAddr[coeffI]] -= mult
        (
            mult.transpose(lower[coeffI]),
            xT[upperAddr[coeffI]]
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockILUC0Precon<Type>::BlockILUC0Precon
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
:
    BlockLduPrecon<Type>(matrix),
    preconDiag_(matrix.diag()),
    preconLower_(matrix.lower()),
    preconUpper_(matrix.upper())
{
    this->calcFactorization();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockILUC0Precon<Type>::~BlockILUC0Precon()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::BlockILUC0Precon<Type>::precondition
(
    Field<Type>& x,
    const Field<Type>& b
) const
{
    if (!this->matrix_.diagonal())
    {
        // Get upper and lower matrix factors
        const CoeffField<Type>& Lower = preconLower_;
        const CoeffField<Type>& Upper = preconUpper_;

        // Execute preconditioning by LU substitution.
        // Note: lower, diag and upper must have same type as required by the
        // algorithm. This is handled by lowest possible promotion
        if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            if (Upper.activeType() == blockCoeffBase::SCALAR)
            {
                LUSubstitute
                (
                    x,
                    preconDiag_.asScalar(),
                    Upper.asScalar(),
                    Lower.asScalar(),
                    b
                );
            }
            else if (Upper.activeType() == blockCoeffBase::LINEAR)
            {
                LUSubstitute
                (
                    x,
                    preconDiag_.asLinear(), // Promotes to linear
                    Upper.asLinear(),
                    Lower.asLinear(),
                    b
                );
            }
            else if (Upper.activeType() == blockCoeffBase::SQUARE)
            {
                LUSubstitute
                (
                    x,
                    preconDiag_.asSquare(), // Promotes to square
                    Upper.asSquare(),
                    Lower.asSquare(),
                    b
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            if (Upper.activeType() == blockCoeffBase::SCALAR)
            {
                LUSubstitute
                (
                    x,
                    preconDiag_.asLinear(),
                    Upper.asLinear(), // Promotes to linear
                    Lower.asLinear(), // Promotes to linear
                    b
                );
            }
            else if (Upper.activeType() == blockCoeffBase::LINEAR)
            {
                LUSubstitute
                (
                    x,
                    preconDiag_.asLinear(),
                    Upper.asLinear(),
                    Lower.asLinear(),
                    b
                );
            }
            else if (Upper.activeType() == blockCoeffBase::SQUARE)
            {
                LUSubstitute
                (
                    x,
                    preconDiag_.asSquare(), // Promotes to square
                    Upper.asSquare(),
                    Lower.asSquare(),
                    b
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::SQUARE)
        {
            if (Upper.activeType() == blockCoeffBase::SCALAR)
            {
                LUSubstitute
                (
                    x,
                    preconDiag_.asSquare(),
                    Upper.asSquare(), // Promotes to square
                    Lower.asSquare(), // Promotes to square
                    b
                );
            }
            else if (Upper.activeType() == blockCoeffBase::LINEAR)
            {
                LUSubstitute
                (
                    x,
                    preconDiag_.asSquare(),
                    Upper.asSquare(), // Promotes to square
                    Lower.asSquare(), // Promotes to square
                    b
                );
            }
            else if (Upper.activeType() == blockCoeffBase::SQUARE)
            {
                LUSubstitute
                (
                    x,
                    preconDiag_.asSquare(),
                    Upper.asSquare(),
                    Lower.asSquare(),
                    b
                );
            }
        }
        else
        {
            FatalErrorIn
            (
                "void BlockILUC0Precon<Type>::precondition\n"
                "(\n"
                "    Field<Type>& x,\n"
                "    const Field<Type>& T\n"
                ") const"
            )   << "Problem with coefficient type morphing."
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn
        (
            "void BlockILUC0Precon<Type>::precondition\n"
            "(\n"
            "    Field<Type>& x,\n"
            "    const Field<Type>& b\n"
            ") const"
        )   << "Unnecessary use of BlockILUC0 preconditioner for diagonal "
            << "matrix. "
            << nl
            << "Use BlockDiagonal preconditioner instead."
            << abort(FatalError);
    }
}


template<class Type>
void Foam::BlockILUC0Precon<Type>::preconditionT
(
    Field<Type>& xT,
    const Field<Type>& bT
) const
{
    if (!this->matrix_.diagonal())
    {
        // Get upper and lower matrix factors
        const CoeffField<Type>& Lower = preconLower_;
        const CoeffField<Type>& Upper = preconUpper_;

        // Execute transpose preconditioning by transpose LU substitution.
        // Note: lower, diag and upper must have same type as required by the
        // algorithm. This is handled by lowest possible promotion
        if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            if (Upper.activeType() == blockCoeffBase::SCALAR)
            {
                LUSubstituteT
                (
                    xT,
                    preconDiag_.asScalar(),
                    Upper.asScalar(),
                    Lower.asScalar(),
                    bT
                );
            }
            else if (Upper.activeType() == blockCoeffBase::LINEAR)
            {
                LUSubstituteT
                (
                    xT,
                    preconDiag_.asLinear(), // Promotes to linear
                    Upper.asLinear(),
                    Lower.asLinear(),
                    bT
                );
            }
            else if (Upper.activeType() == blockCoeffBase::SQUARE)
            {
                LUSubstituteT
                (
                    xT,
                    preconDiag_.asSquare(), // Promotes to square
                    Upper.asSquare(),
                    Lower.asSquare(),
                    bT
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            if (Upper.activeType() == blockCoeffBase::SCALAR)
            {
                LUSubstituteT
                (
                    xT,
                    preconDiag_.asLinear(),
                    Upper.asLinear(), // Promotes to linear
                    Lower.asLinear(), // Promotes to linear
                    bT
                );
            }
            else if (Upper.activeType() == blockCoeffBase::LINEAR)
            {
                LUSubstituteT
                (
                    xT,
                    preconDiag_.asLinear(),
                    Upper.asLinear(),
                    Lower.asLinear(),
                    bT
                );
            }
            else if (Upper.activeType() == blockCoeffBase::SQUARE)
            {
                LUSubstituteT
                (
                    xT,
                    preconDiag_.asSquare(), // Promotes to square
                    Upper.asSquare(),
                    Lower.asSquare(),
                    bT
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::SQUARE)
        {
            if (Upper.activeType() == blockCoeffBase::SCALAR)
            {
                LUSubstituteT
                (
                    xT,
                    preconDiag_.asSquare(),
                    Upper.asSquare(), // Promotes to square
                    Lower.asSquare(), // Promotes to square
                    bT
                );
            }
            else if (Upper.activeType() == blockCoeffBase::LINEAR)
            {
                LUSubstituteT
                (
                    xT,
                    preconDiag_.asSquare(),
                    Upper.asSquare(), // Promotes to square
                    Lower.asSquare(), // Promotes to square
                    bT
                );
            }
            else if (Upper.activeType() == blockCoeffBase::SQUARE)
            {
                LUSubstituteT
                (
                    xT,
                    preconDiag_.asSquare(),
                    Upper.asSquare(),
                    Lower.asSquare(),
                    bT
                );
            }
        }
        else
        {
            FatalErrorIn
            (
                "void BlockILUC0Precon<Type>::preconditionT\n"
                "(\n"
                "    Field<Type>& x,\n"
                "    const Field<Type>& T\n"
                ") const"
            )   << "Problem with coefficient type morphing."
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn
        (
            "void BlockILUC0Precon<Type>::preconditionT\n"
            "(\n"
            "    Field<Type>& x,\n"
            "    const Field<Type>& b\n"
            ") const"
        )   << "Unnecessary use of BlockILUC0 preconditioner for diagonal "
            << "matrix."
            << nl
            << "Use BlockDiagonal preconditioner instead."
            << abort(FatalError);
    }
}


template<class Type>
void Foam::BlockILUC0Precon<Type>::initMatrix()
{
    preconDiag_ = this->matrix_.diag();
    preconLower_ = this->matrix_.lower();
    preconUpper_ = this->matrix_.upper();

    this->calcFactorization();
}


// ************************************************************************* //
