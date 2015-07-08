/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "BlockILUCpPrecon.H"
#include "error.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class LDUType>
void Foam::BlockILUCpPrecon<Type>::calcFactorization
(
    Field<LDUType>& preconD,
    Field<LDUType>& extUpper,
    Field<LDUType>& extLower,
    Field<LDUType>& zDiag,
    Field<LDUType>& z,
    Field<LDUType>& w
)
{
    if (!this->matrix_.diagonal())
    {
        // Create multiplication function object
        typename BlockCoeff<LDUType>::multiply mult;

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
        LDUType* __restrict__ diagPtr = preconD.begin();
        LDUType* __restrict__ upperPtr = extUpper.begin();
        LDUType* __restrict__ lowerPtr = extLower.begin();

        // Get access to working fields
        LDUType* __restrict__ zPtr = z.begin();
        LDUType* __restrict__ wPtr = w.begin();
        LDUType& zDiagI = zDiag[0];

        // Get number of rows
        const label nRows = preconD.size();

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
            zDiagI = diagPtr[rowI];

            // Initialize temporary working row and column fields
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
                // WARNING: Not sure about order of multiplication.
                // Check it. VV, 3/Jul/2015
                zDiagI -= mult
                (
                    lowerPtr[losortCoeff],
                    upperPtr[losortCoeff]
                );

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
                    // WARNING: Not sure about order of multiplication.
                    // Check it. VV, 3/Jul/2015
                    zPtr[uPtr[faceI]] -= mult
                    (
                        lowerPtr[losortCoeff],
                        upperPtr[faceI]
                    );

                    // WARNING: Not sure about order of multiplication.
                    // Check it. VV, 3/Jul/2015
                    wPtr[uPtr[faceI]] -= mult
                    (
                        upperPtr[losortCoeff],
                        lowerPtr[faceI]
                    );
                }
            }

            // Update diagonal entry, inverting it for future use
            LDUType& diagRowI = diagPtr[rowI];
            diagRowI = mult.inverse(zDiagI);

            // Index for updating L and U
            register label zwIndex;

            // Update upper and lower coeffs
            for (register label faceI = fStart; faceI < fEnd; ++faceI)
            {
                // Get index for current face
                zwIndex = uPtr[faceI];

                // Update L and U decomposition for this row (column)
                upperPtr[faceI] = zPtr[zwIndex];

                // WARNING: Not sure about order of multiplication.
                // Check it. VV, 3/Jul/2015
                lowerPtr[faceI] = mult
                (
                    diagRowI,
                    wPtr[zwIndex]
                );
            }

            // Reset temporary working fields
            zDiagI = pTraits<LDUType>::zero;

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
            "void BlockILUCpPrecon<Type>::calcFactorization\n"
            "(\n"
            "    Field<LDUType>& preconD,\n"
            "    Field<LDUType>& extUpper,\n"
            "    Field<LDUType>& extLower,\n"
            "    field<LDUType>& zDiag\n,"
            "    Field<LDUType>& z,\n"
            "    Field<LDUType>& w,\n"
            ") const"
        )   << "Unnecessary use of BlockILUCp preconditioner for diagonal "
            << "matrix."
            << nl
            << "Use BlockDiagonal preconditioner instead."
            << abort(FatalError);
    }
}


template<class Type>
template<class LDUType>
void Foam::BlockILUCpPrecon<Type>::LUSubstitute
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
    const extendedLduAddressing& addr = extBlockMatrix_.extendedLduAddr();
    const unallocLabelList& upperAddr = addr.extendedUpperAddr();
    const unallocLabelList& lowerAddr = addr.extendedLowerAddr();
    const unallocLabelList& losortAddr = addr.extendedLosortAddr();

    // Solve Lz = b with forward substitution in block form. lower is chosen
    // to be unit triangular. z does not need to be stored

    // Initialize x field
    x = b;

    register label losortCoeffI;
    register label rowI;

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
void Foam::BlockILUCpPrecon<Type>::LUSubstituteT
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
    const extendedLduAddressing& addr = extBlockMatrix_.extendedLduAddr();
    const unallocLabelList& upperAddr = addr.extendedUpperAddr();
    const unallocLabelList& lowerAddr = addr.extendedLowerAddr();
    const unallocLabelList& losortAddr = addr.extendedLosortAddr();

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

    register label losortCoeffI;
    register label rowI;

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
Foam::BlockILUCpPrecon<Type>::BlockILUCpPrecon
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
:
    BlockLduPrecon<Type>(matrix),
    preconDiag_(matrix.diag()),
    p_(readLabel(dict.lookup("fillInLevel"))),
    extBlockMatrix_
    (
        matrix,
        p_,
        matrix.mesh().thisDb().time().template lookupObject<polyMesh>
        (
            polyMesh::defaultRegion
        )
    ),
    zDiag_(1),
    z_(preconDiag_.size()),
    w_(preconDiag_.size())
{
    // Get upper and lower matrix factors
    CoeffField<Type>& Lower = extBlockMatrix_.extendedLower();
    CoeffField<Type>& Upper = extBlockMatrix_.extendedUpper();

    // Calculate factorization
    // Note: lower, diag and upper must have same type as required by the
    // algorithm. This is handled by lowest possible promotion
    if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
    {
        if (Upper.activeType() == blockCoeffBase::SCALAR)
        {
            calcFactorization
            (
                preconDiag_.asScalar(),
                Upper.asScalar(),
                Lower.asScalar(),
                zDiag_.asScalar(),
                z_.asScalar(),
                w_.asScalar()
            );
        }
        else if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
            calcFactorization
            (
                preconDiag_.asLinear(), // Promotes to linear
                Upper.asLinear(),
                Lower.asLinear(),
                zDiag_.asLinear(),
                z_.asLinear(),
                w_.asLinear()
            );
        }
        else if (Upper.activeType() == blockCoeffBase::SQUARE)
        {
            calcFactorization
            (
                preconDiag_.asSquare(), // Promotes to square
                Upper.asSquare(),
                Lower.asSquare(),
                zDiag_.asSquare(),
                z_.asSquare(),
                w_.asSquare()
            );
        }
    }
    else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
    {
        if (Upper.activeType() == blockCoeffBase::SCALAR)
        {
            calcFactorization
            (
                preconDiag_.asLinear(),
                Upper.asLinear(), // Promotes to linear
                Lower.asLinear(), // Promotes to linear
                zDiag_.asLinear(), // Promotes to linear
                z_.asLinear(), // Promotes to linear
                w_.asLinear() // Promotes to linear
            );
        }
        else if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
            calcFactorization
            (
                preconDiag_.asLinear(),
                Upper.asLinear(),
                Lower.asLinear(),
                zDiag_.asLinear(),
                z_.asLinear(),
                w_.asLinear()
            );
        }
        else if (Upper.activeType() == blockCoeffBase::SQUARE)
        {
            calcFactorization
            (
                preconDiag_.asSquare(), // Promotes to square
                Upper.asSquare(),
                Lower.asSquare(),
                zDiag_.asSquare(),
                z_.asSquare(),
                w_.asSquare()
            );
        }
    }
    else if (preconDiag_.activeType() == blockCoeffBase::SQUARE)
    {
        if (Upper.activeType() == blockCoeffBase::SCALAR)
        {
            calcFactorization
            (
                preconDiag_.asSquare(),
                Upper.asSquare(), // Promotes to square
                Lower.asSquare(), // Promotes to square
                zDiag_.asSquare(), // Promotes to square
                z_.asSquare(), // Promotes to square
                w_.asSquare() // Promotes to square
            );
        }
        else if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
            calcFactorization
            (
                preconDiag_.asSquare(),
                Upper.asSquare(), // Promotes to square
                Lower.asSquare(), // Promotes to square
                zDiag_.asSquare(), // Promotes to square
                z_.asSquare(), // Promotes to square
                w_.asSquare() // Promotes to square
            );
        }
        else if (Upper.activeType() == blockCoeffBase::SQUARE)
        {
            calcFactorization
            (
                preconDiag_.asSquare(),
                Upper.asSquare(),
                Lower.asSquare(),
                zDiag_.asSquare(),
                z_.asSquare(),
                w_.asSquare()
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockILUCpPrecon<Type>::~BlockILUCpPrecon()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::BlockILUCpPrecon<Type>::precondition
(
    Field<Type>& x,
    const Field<Type>& b
) const
{
    if (!this->matrix_.diagonal())
    {
        // Get upper and lower matrix factors
        const CoeffField<Type>& Lower = extBlockMatrix_.extendedLower();
        const CoeffField<Type>& Upper = extBlockMatrix_.extendedUpper();

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
                "void BlockILUCpPrecon<Type>::precondition\n"
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
            "void BlockILUCpPrecon<Type>::precondition\n"
            "(\n"
            "    Field<Type>& x,\n"
            "    const Field<Type>& b\n"
            ") const"
        )   << "Unnecessary use of BlockILUCp preconditioner for diagonal "
            << "matrix. "
            << nl
            << "Use BlockDiagonal preconditioner instead."
            << abort(FatalError);
    }
}


template<class Type>
void Foam::BlockILUCpPrecon<Type>::preconditionT
(
    Field<Type>& xT,
    const Field<Type>& bT
) const
{
    if (!this->matrix_.diagonal())
    {
        // Get upper and lower matrix factors
        const CoeffField<Type>& Lower = extBlockMatrix_.extendedLower();
        const CoeffField<Type>& Upper = extBlockMatrix_.extendedUpper();

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
                "void BlockILUCpPrecon<Type>::preconditionT\n"
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
            "void BlockILUCpPrecon<Type>::preconditionT\n"
            "(\n"
            "    Field<Type>& x,\n"
            "    const Field<Type>& b\n"
            ") const"
        )   << "Unnecessary use of BlockILUCp preconditioner for diagonal "
            << "matrix."
            << nl
            << "Use BlockDiagonal preconditioner instead."
            << abort(FatalError);
    }
}


// ************************************************************************* //
