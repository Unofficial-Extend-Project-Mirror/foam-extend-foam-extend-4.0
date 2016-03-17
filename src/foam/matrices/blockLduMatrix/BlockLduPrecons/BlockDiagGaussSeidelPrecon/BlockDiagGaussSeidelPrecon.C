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

Description
    Diagonally-corrected block Gauss-Seidel sweep as a preconditioner.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "BlockDiagGaussSeidelPrecon.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::BlockDiagGaussSeidelPrecon<Type>::calcInvDiag()
{
    typedef CoeffField<Type> TypeCoeffField;
    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
    typedef typename TypeCoeffField::scalarType scalarType;

    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::linearType linearType;

    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    const unallocLabelList& l = this->matrix_.lduAddr().lowerAddr();
    const unallocLabelList& u = this->matrix_.lduAddr().upperAddr();

    // Get reference to diagonal
    const TypeCoeffField& d = this->matrix_.diag();

    TypeCoeffField sumMagOffDiag(d.size());

    // Sum up off-diagonal part of the matrix

    if (this->matrix_.symmetric())
    {
        // Symmetric matrix

        const TypeCoeffField& Upper = this->matrix_.upper();

        if (Upper.activeType() == blockCoeffBase::SQUARE)
        {
            // magOffDiag is always scalar or linear type
            linearTypeField& activeDiag = sumMagOffDiag.asLinear();

            // Use lower as transpose of upper
            linearTypeField activeLower
            (
                Upper.size(),
                pTraits<linearType>::zero
            );

            sumMagToDiag(activeLower, Upper.asSquare().T()());

            linearTypeField activeUpper
            (
                Upper.size(),
                pTraits<linearType>::zero
            );

            sumMagToDiag(activeUpper, Upper.asSquare());

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeDiag[l[coeffI]] += activeLower[coeffI];
                activeDiag[u[coeffI]] += activeUpper[coeffI];
            }
        }
        else if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
            // magOffDiag is always scalar or linear type
            linearTypeField& activeDiag = sumMagOffDiag.asLinear();

            // Lower is identical to upper
            const linearTypeField& activeUpper = Upper.asLinear();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeDiag[l[coeffI]] += activeUpper[coeffI];
                activeDiag[u[coeffI]] += activeUpper[coeffI];
            }
        }
        else if (Upper.activeType() == blockCoeffBase::SCALAR)
        {
            // magOffDiag is always scalar or linear type
            scalarTypeField& activeDiag = sumMagOffDiag.asScalar();

            // Lower is identical to upper
            const scalarTypeField& activeUpper = Upper.asScalar();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeDiag[l[coeffI]] += activeUpper[coeffI];
                activeDiag[u[coeffI]] += activeUpper[coeffI];
            }
        }
    }
    else if (this->matrix_.asymmetric())
    {
        // Assymetric matrix

        const TypeCoeffField& Lower = this->matrix_.lower();
        const TypeCoeffField& Upper = this->matrix_.upper();

        if (Upper.activeType() == blockCoeffBase::SQUARE)
        {
            // magOffDiag is always scalar or linear type
            linearTypeField& activeDiag = sumMagOffDiag.asLinear();

            linearTypeField activeLower
            (
                Lower.size(),
                pTraits<linearType>::zero
            );

            sumMagToDiag(activeLower, Lower.asSquare());

            linearTypeField activeUpper
            (
                Upper.size(),
                pTraits<linearType>::zero
            );

            sumMagToDiag(activeUpper, Upper.asSquare());

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeDiag[l[coeffI]] += activeLower[coeffI];
                activeDiag[u[coeffI]] += activeUpper[coeffI];
            }
        }
        else if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
            // magOffDiag is always scalar or linear type
            linearTypeField& activeDiag = sumMagOffDiag.asLinear();

            const linearTypeField& activeLower = Lower.asLinear();
            const linearTypeField& activeUpper = Upper.asLinear();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeDiag[l[coeffI]] += activeLower[coeffI];
                activeDiag[u[coeffI]] += activeUpper[coeffI];
            }
        }
        else if (Upper.activeType() == blockCoeffBase::SCALAR)
        {
            scalarTypeField& activeDiag = sumMagOffDiag.asScalar();

            const scalarTypeField& activeLower = Lower.asScalar();
            const scalarTypeField& activeUpper = Upper.asScalar();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeDiag[l[coeffI]] += activeLower[coeffI];
                activeDiag[u[coeffI]] += activeUpper[coeffI];
            }
        }
    }


    // Compare sum off-diag with diagonal and take larger for the inverse
    // The difference between the two is stored in LUDiag_ for correction

    if (d.activeType() == blockCoeffBase::SCALAR)
    {
        if (sumMagOffDiag.activeType() == blockCoeffBase::SCALAR)
        {
            const scalarTypeField& activeDiag = d.asScalar();
            const scalarTypeField& activeCorrDiag = sumMagOffDiag.asScalar();

            scalarTypeField signActiveDiag = sign(activeDiag);
            scalarTypeField magActiveDiag = mag(activeDiag);

            scalarTypeField magNewDiag =
                Foam::max(magActiveDiag, activeCorrDiag);

            // Calculate inverse with corrected diagonal
            invDiag_.asScalar() = signActiveDiag/magNewDiag;

            // Store correction into LUDiag_
            LUDiag_.asScalar() = signActiveDiag*
                Foam::max
                (
                    magNewDiag - magActiveDiag,
                    pTraits<scalarType>::zero
                );
        }
        else if (sumMagOffDiag.activeType() == blockCoeffBase::LINEAR)
        {
            const scalarTypeField& activeDiag = d.asScalar();
            const linearTypeField& activeCorrDiag = sumMagOffDiag.asLinear();

            scalarTypeField signActiveDiag = sign(activeDiag);
            linearTypeField magActiveDiag =
                mag(activeDiag)*pTraits<linearType>::one;

            linearTypeField magNewDiag =
                Foam::max(magActiveDiag, activeCorrDiag);

            invDiag_.asLinear() = signActiveDiag*
                cmptDivide
                (
                    pTraits<linearType>::one,
                    magNewDiag
                );

            // Store correction into LUDiag_
            LUDiag_.asLinear() = signActiveDiag*
                max
                (
                    magNewDiag - magActiveDiag,
                    pTraits<linearType>::zero
                );
        }
    }
    else if (d.activeType() == blockCoeffBase::LINEAR)
    {
        const linearTypeField& activeDiag = d.asLinear();
        const linearTypeField& activeCorrDiag = sumMagOffDiag.asLinear();

        linearTypeField signActiveDiag = cmptSign(activeDiag);
        linearTypeField magActiveDiag = cmptMag(activeDiag);

        linearTypeField magNewDiag =
            Foam::max(magActiveDiag, activeCorrDiag);

        invDiag_.asLinear() =
            cmptDivide
            (
                signActiveDiag,
                magNewDiag
            );

        // Store correction into LUDiag_
        LUDiag_.asLinear() =
            cmptMultiply
            (
                signActiveDiag,
                max
                (
                    magNewDiag - magActiveDiag,
                    pTraits<linearType>::zero
                )
            );
    }
    else if (d.activeType() == blockCoeffBase::SQUARE)
    {
        // For square diagonal invert diagonal only and store the rest
        // info LUDiag coefficient.  This avoids full inverse of the
        // diagonal matrix.  HJ, 20/Aug/2015

        // Get reference to active diag
        const squareTypeField& activeDiag = d.asSquare();

        // Get reference to LU: remove diagonal from active diag
        squareTypeField& luDiag = LUDiag_.asSquare();

        linearTypeField diagDiag(activeDiag.size());

        // Take out the diagonal from the diag as a linear type
        contractLinear(diagDiag, activeDiag);

        // Expand diagonal only to full square type and store into luDiag
        expandLinear(luDiag, diagDiag);

        // Keep only off-diagonal in luDiag.
        // Note change of sign to avoid multiplication with -1 when moving
        // to the other side.  HJ, 20/Aug/2015
        luDiag -= activeDiag;

        // The diagonal is now in diagDiag

        const linearTypeField& activeCorrDiag = sumMagOffDiag.asLinear();
        linearTypeField signActiveDiag = cmptSign(diagDiag);
        linearTypeField magActiveDiag = cmptMag(diagDiag);

        linearTypeField magNewDiag =
            Foam::max(magActiveDiag, activeCorrDiag);

        invDiag_.asLinear() =
            cmptDivide
            (
                signActiveDiag,
                magNewDiag
            );

        // Inverse is the inverse of diagonal only
        linearTypeField diagDiagCorr =
            cmptMultiply
            (
                signActiveDiag,
                Foam::max
                (
                    magNewDiag - magActiveDiag,
                    pTraits<linearType>::zero
                )
            );

        // Add correction to luDiag
        squareTypeField corr(d.size());

        expandLinear(corr, diagDiagCorr);

        luDiag += corr;
    }
    else
    {
        FatalErrorIn
        (
            "void BlockDiagGaussSeidelPrecon<Type>::calcInvDiag()"
        )   << "Problem with coefficient type morphing."
            << abort(FatalError);
    }
}


// Block sweep, symmetric matrix
template<class Type>
template<class DiagType, class ULType>
void Foam::BlockDiagGaussSeidelPrecon<Type>::BlockSweep
(
    Field<Type>& x,
    const Field<DiagType>& dD,
    const Field<ULType>& upper,
    const Field<Type>& b
) const
{
    const unallocLabelList& u = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& ownStart =
        this->matrix_.lduAddr().ownerStartAddr();

    const label nRows = ownStart.size() - 1;

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    // Klas Jareteg: 2013-02-10:
    // Must transfer data between the different CPUs. Notes on the Jacobi
    // iteration style can be seen in DiagGaussSeidelSolver.C

    bPrime_ = b;

    this->matrix_.initInterfaces
    (
        this->matrix_.coupleUpper(),
        bPrime_,
        x,
        true             // switch to lhs of system
    );

    this->matrix_.updateInterfaces
    (
        this->matrix_.coupleUpper(),
        bPrime_,
        x,
        true             // switch to lhs of system
    );

    register label fStart, fEnd, curCoeff;

    // Forward sweep
    for (register label rowI = 0; rowI < nRows; rowI++)
    {
        Type& curX = x[rowI];

        // Grab the accumulated neighbour side
        curX = bPrime_[rowI];

        // Start and end of this row
        fStart = ownStart[rowI];
        fEnd = ownStart[rowI + 1];

        // Accumulate the owner product side
        for (curCoeff = fStart; curCoeff < fEnd; curCoeff++)
        {
            curX -= mult(upper[curCoeff], x[u[curCoeff]]);
        }

        // Finish current x
        curX = mult(dD[rowI], curX);

        // Distribute the neighbour side using current x
        for (curCoeff = fStart; curCoeff < fEnd; curCoeff++)
        {
            // lower = upper transposed
            bPrime_[u[curCoeff]] -=
                mult(mult.transpose(upper[curCoeff]), curX);
        }
    }

    // Reverse sweep
    for (register label rowI = nRows - 1; rowI >= 0; rowI--)
    {
        Type& curX = x[rowI];

        // Grab the accumulated neighbour side
        curX = bPrime_[rowI];

        // Start and end of this row
        fStart = ownStart[rowI];
        fEnd = ownStart[rowI + 1];

        // Accumulate the owner product side
        for (curCoeff = fStart; curCoeff < fEnd; curCoeff++)
        {
            curX -= mult(upper[curCoeff], x[u[curCoeff]]);
        }

        // Finish current x
        curX = mult(dD[rowI], curX);

        // No need to update bPrime on reverse sweep. VV, 10/Sep/2015.
    }
}


// Block sweep, asymmetric matrix
template<class Type>
template<class DiagType, class ULType>
void Foam::BlockDiagGaussSeidelPrecon<Type>::BlockSweep
(
    Field<Type>& x,
    const Field<DiagType>& dD,
    const Field<ULType>& lower,
    const Field<ULType>& upper,
    const Field<Type>& b
) const
{
    const unallocLabelList& u = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& ownStart =
        this->matrix_.lduAddr().ownerStartAddr();

    const label nRows = ownStart.size() - 1;

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    // Klas Jareteg: 2013-02-10:
    // Must transfer data between the different CPUs. Notes on the Jacobi
    // iteration style can be seen in DiagGaussSeidelSolver.C

    bPrime_ = b;

    this->matrix_.initInterfaces
    (
        this->matrix_.coupleUpper(),
        bPrime_,
        x,
        true             // switch to lhs of system
    );

    this->matrix_.updateInterfaces
    (
        this->matrix_.coupleUpper(),
        bPrime_,
        x,
        true             // switch to lhs of system
    );

    register label fStart, fEnd, curCoeff;

    // Forward sweep
    for (register label rowI = 0; rowI < nRows; rowI++)
    {
        Type& curX = x[rowI];

        // Grab the accumulated neighbour side
        curX = bPrime_[rowI];

        // Start and end of this row
        fStart = ownStart[rowI];
        fEnd = ownStart[rowI + 1];

        // Accumulate the owner product side
        for (curCoeff = fStart; curCoeff < fEnd; curCoeff++)
        {
            curX -= mult(upper[curCoeff], x[u[curCoeff]]);
        }

        // Finish current x
        curX = mult(dD[rowI], curX);

        // Distribute the neighbour side using current x
        for (curCoeff = fStart; curCoeff < fEnd; curCoeff++)
        {
            bPrime_[u[curCoeff]] -= mult(lower[curCoeff], curX);
        }
    }

    // Reverse sweep
    for (register label rowI = nRows - 1; rowI >= 0; rowI--)
    {
        Type& curX = x[rowI];

        // Grab the accumulated neighbour side
        curX = bPrime_[rowI];

        // Start and end of this row
        fStart = ownStart[rowI];
        fEnd = ownStart[rowI + 1];

        // Accumulate the owner product side
        for (curCoeff = fStart; curCoeff < fEnd; curCoeff++)
        {
            curX -= mult(upper[curCoeff], x[u[curCoeff]]);
        }

        // Finish current x
        curX = mult(dD[rowI], curX);

        // No need to update bPrime on reverse sweep. VV, 10/Sep/2015.
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockDiagGaussSeidelPrecon<Type>::BlockDiagGaussSeidelPrecon
(
    const BlockLduMatrix<Type>& matrix
)
:
    BlockLduPrecon<Type>(matrix),
    invDiag_(matrix.lduAddr().size()),
    LUDiag_(matrix.lduAddr().size()),
    bPlusLU_(),
    bPrime_(matrix.lduAddr().size()),
    nSweeps_(1)
{
    calcInvDiag();
}


template<class Type>
Foam::BlockDiagGaussSeidelPrecon<Type>::BlockDiagGaussSeidelPrecon
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
:
    BlockLduPrecon<Type>(matrix),
    invDiag_(matrix.lduAddr().size()),
    LUDiag_(matrix.lduAddr().size()),
    bPlusLU_(),
    bPrime_(matrix.lduAddr().size()),
    nSweeps_(readLabel(dict.lookup("nSweeps")))
{
    calcInvDiag();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::BlockDiagGaussSeidelPrecon<Type>::precondition
(
    Field<Type>& x,
    const Field<Type>& b
) const
{
    typedef CoeffField<Type> TypeCoeffField;

    // Prepare correction
    if (bPlusLU_.empty())
    {
        bPlusLU_.setSize(b.size(), pTraits<Type>::zero);
    }

    if (this->matrix_.diagonal())
    {
        multiply(x, invDiag_, b);
    }
    else if (this->matrix_.symmetric())
    {
        const TypeCoeffField& DiagCoeff = this->matrix_.diag();
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        // Note
        // Gauss-Seidel loops need to be executed in the specific
        // order with direct access to the coefficients which can be
        // of morphed type.  Under normal circumstances, the
        // operations are not inter-leaved and the decision can be
        // made at the beginning of the loop.  Here, the order needs
        // to be enforced without the per-element if-condition, which
        // makes for ugly code.  HJ, 19/May/2005

        // Note: Assuming lower and upper triangle have the same active type

        if (DiagCoeff.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, x);
                    bPlusLU_ += b;

                    BlockSweep
                    (
                        x,
                        invDiag_.asScalar(),
                        UpperCoeff.asScalar(),
                        bPlusLU_
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, x);
                    bPlusLU_ += b;

                    BlockSweep
                    (
                        x,
                        invDiag_.asScalar(),
                        UpperCoeff.asLinear(),
                        bPlusLU_
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, x);
                    bPlusLU_ += b;

                    BlockSweep
                    (
                        x,
                        invDiag_.asScalar(),
                        UpperCoeff.asSquare(),
                        bPlusLU_
                    );
                }
            }
        }
        else if (DiagCoeff.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, x);
                    bPlusLU_ += b;

                    BlockSweep
                    (
                        x,
                        invDiag_.asLinear(),
                        UpperCoeff.asScalar(),
                        bPlusLU_
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, x);
                    bPlusLU_ += b;

                    BlockSweep
                    (
                        x,
                        invDiag_.asLinear(),
                        UpperCoeff.asLinear(),
                        bPlusLU_
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, x);
                    bPlusLU_ += b;

                    BlockSweep
                    (
                        x,
                        invDiag_.asLinear(),
                        UpperCoeff.asSquare(),
                        bPlusLU_
                    );
                }
            }
        }
        else if (DiagCoeff.activeType() == blockCoeffBase::SQUARE)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, x);
                    bPlusLU_ += b;

                    // Note linear diag inversed due to decoupling
                    BlockSweep
                    (
                        x,
                        invDiag_.asLinear(),
                        UpperCoeff.asScalar(),
                        bPlusLU_
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, x);
                    bPlusLU_ += b;

                    // Note linear diag inversed due to decoupling
                    BlockSweep
                    (
                        x,
                        invDiag_.asLinear(),
                        UpperCoeff.asLinear(),
                        bPlusLU_
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, x);
                    bPlusLU_ += b;

                    // Note linear diag inversed due to decoupling
                    BlockSweep
                    (
                        x,
                        invDiag_.asLinear(),
                        UpperCoeff.asSquare(),
                        bPlusLU_
                    );
                }
            }
        }
        else
        {
            FatalErrorIn
            (
                "void BlockGaussSeidelPrecon<Type>::precondition\n"
                "(\n"
                "    Field<Type>& x,\n"
                "    const Field<Type>& b\n"
                ") const"
             )  << "Problem with coefficient type morphing."
                << abort(FatalError);
        }
    }
    else if (this->matrix_.asymmetric())
    {
        const TypeCoeffField& DiagCoeff = this->matrix_.diag();
        const TypeCoeffField& LowerCoeff = this->matrix_.lower();
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        // Note
        // Gauss-Seidel loops need to be executed in the specific
        // order with direct access to the coefficients which can be
        // of morphed type.  Under normal circumstances, the
        // operations are not inter-leaved and the decision can be
        // made at the beginning of the loop.  Here, the order needs
        // to be enforced without the per-element if-condition, which
        // makes for ugly code.  HJ, 19/May/2005

        //Note: Assuming lower and upper triangle have the same active type

        if (DiagCoeff.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, x);
                    bPlusLU_ += b;

                    BlockSweep
                    (
                        x,
                        invDiag_.asScalar(),
                        LowerCoeff.asScalar(),
                        UpperCoeff.asScalar(),
                        bPlusLU_
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, x);
                    bPlusLU_ += b;

                    BlockSweep
                    (
                        x,
                        invDiag_.asScalar(),
                        LowerCoeff.asLinear(),
                        UpperCoeff.asLinear(),
                        bPlusLU_
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, x);
                    bPlusLU_ += b;

                    BlockSweep
                    (
                        x,
                        invDiag_.asScalar(),
                        LowerCoeff.asSquare(),
                        UpperCoeff.asSquare(),
                        bPlusLU_
                    );
                }
            }
        }
        else if (DiagCoeff.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, x);
                    bPlusLU_ += b;

                    BlockSweep
                    (
                        x,
                        invDiag_.asLinear(),
                        LowerCoeff.asScalar(),
                        UpperCoeff.asScalar(),
                        bPlusLU_
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, x);
                    bPlusLU_ += b;

                    BlockSweep
                    (
                        x,
                        invDiag_.asLinear(),
                        LowerCoeff.asLinear(),
                        UpperCoeff.asLinear(),
                        bPlusLU_
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, x);
                    bPlusLU_ += b;

                    BlockSweep
                    (
                        x,
                        invDiag_.asLinear(),
                        LowerCoeff.asSquare(),
                        UpperCoeff.asSquare(),
                        bPlusLU_
                    );
                }
            }
        }
        else if (DiagCoeff.activeType() == blockCoeffBase::SQUARE)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, x);
                    bPlusLU_ += b;

                    // Note linear diag inversed due to decoupling
                    BlockSweep
                    (
                        x,
                        invDiag_.asLinear(),
                        LowerCoeff.asScalar(),
                        UpperCoeff.asScalar(),
                        bPlusLU_
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, x);
                    bPlusLU_ += b;

                    // Note linear diag inversed due to decoupling
                    BlockSweep
                    (
                        x,
                        invDiag_.asLinear(),
                        LowerCoeff.asLinear(),
                        UpperCoeff.asLinear(),
                        bPlusLU_
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, x);
                    bPlusLU_ += b;

                    // Note linear diag inversed due to decoupling
                    BlockSweep
                    (
                        x,
                        invDiag_.asLinear(),
                        LowerCoeff.asSquare(),
                        UpperCoeff.asSquare(),
                        bPlusLU_
                    );
                }
            }
        }
        else
        {
            FatalErrorIn
            (
                "void BlockDiagGaussSeidelPrecon<Type>::precondition\n"
                "(\n"
                "    Field<Type>& x,\n"
                "    const Field<Type>& b\n"
                ") const"
             )  << "Problem with coefficient type morphing."
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn
        (
            "void BlockDiagGaussSeidelPrecon<Type>::precondition\n"
            "(\n"
            "    Field<Type>& x,\n"
            "    const Field<Type>& b\n"
            ") const"
         )  << "cannot solve incomplete matrix, no diagonal"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::BlockDiagGaussSeidelPrecon<Type>::preconditionT
(
    Field<Type>& xT,
    const Field<Type>& bT
) const
{
    typedef CoeffField<Type> TypeCoeffField;

    if (bPlusLU_.empty())
    {
        bPlusLU_.setSize(bT.size(), pTraits<Type>::zero);
    }

    if (this->matrix_.diagonal())
    {
        multiply(xT, invDiag_, bT);
    }
    else if (this->matrix_.symmetric() || this->matrix_.asymmetric())
    {
        const TypeCoeffField& DiagCoeff = this->matrix_.lower();
        const TypeCoeffField& LowerCoeff = this->matrix_.lower();
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        // Note
        // Gauss-Seidel loops need to be executed in the specific
        // order with direct access to the coefficients which can be
        // of morphed type.  Under normal circumstances, the
        // operations are not inter-leaved and the decision can be
        // made at the beginning of the loop.  Here, the order needs
        // to be enforced without the per-element if-condition, which
        // makes for ugly code.  HJ, 19/May/2005

        // Note: Assuming lower and upper triangle have the same active type

        if (DiagCoeff.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, xT);
                    bPlusLU_ += bT;

                    // Transpose multiplication - swap lower and upper coeff
                    BlockSweep
                    (
                        xT,
                        invDiag_.asScalar(),
                        UpperCoeff.asScalar(),
                        LowerCoeff.asScalar(),
                        bPlusLU_
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, xT);
                    bPlusLU_ += bT;

                    // Transpose multiplication - swap lower and upper coeff
                    BlockSweep
                    (
                        xT,
                        invDiag_.asScalar(),
                        UpperCoeff.asLinear(),
                        LowerCoeff.asLinear(),
                        bPlusLU_
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, xT);
                    bPlusLU_ += bT;

                    // Transpose multiplication - swap lower and upper coeff
                    BlockSweep
                    (
                        xT,
                        invDiag_.asScalar(),
                        UpperCoeff.asSquare(),
                        LowerCoeff.asSquare(),
                        bPlusLU_
                    );
                }
            }
        }
        else if (DiagCoeff.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, xT);
                    bPlusLU_ += bT;

                    // Transpose multiplication - swap lower and upper coeff
                    BlockSweep
                    (
                        xT,
                        invDiag_.asLinear(),
                        UpperCoeff.asScalar(),
                        LowerCoeff.asScalar(),
                        bPlusLU_
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, xT);
                    bPlusLU_ += bT;

                    // Transpose multiplication - swap lower and upper coeff
                    BlockSweep
                    (
                        xT,
                        invDiag_.asLinear(),
                        UpperCoeff.asLinear(),
                        LowerCoeff.asLinear(),
                        bPlusLU_
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, xT);
                    bPlusLU_ += bT;

                    // Transpose multiplication - swap lower and upper coeff
                    BlockSweep
                    (
                        xT,
                        invDiag_.asLinear(),
                        UpperCoeff.asSquare(),
                        LowerCoeff.asSquare(),
                        bPlusLU_
                    );
                }
            }
        }
        else if (DiagCoeff.activeType() == blockCoeffBase::SQUARE)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, xT);
                    bPlusLU_ += bT;

                    // Transpose multiplication - swap lower and upper coeff
                    // Note linear diag inversed due to decoupling
                    BlockSweep
                    (
                        xT,
                        invDiag_.asLinear(),
                        UpperCoeff.asScalar(),
                        LowerCoeff.asScalar(),
                        bPlusLU_
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, xT);
                    bPlusLU_ += bT;

                    // Transpose multiplication - swap lower and upper coeff
                    // Note linear diag inversed due to decoupling
                    BlockSweep
                    (
                        xT,
                        invDiag_.asLinear(),
                        UpperCoeff.asLinear(),
                        LowerCoeff.asLinear(),
                        bPlusLU_
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Add diag coupling to b
                    multiply(bPlusLU_, LUDiag_, xT);
                    bPlusLU_ += bT;

                    // Transpose multiplication - swap lower and upper coeff
                    // Note linear diag inversed due to decoupling
                    BlockSweep
                    (
                        xT,
                        invDiag_.asLinear(),
                        UpperCoeff.asSquare(),
                        LowerCoeff.asSquare(),
                        bPlusLU_
                    );
                }
            }
        }
        else
        {
            FatalErrorIn
            (
                "void BlockDiagGaussSeidelPrecon<Type>::preconditionT\n"
                "(\n"
                "    Field<Type>& xT,\n"
                "    const Field<Type>& bT\n"
                ") const"
             )  << "Problem with coefficient type morphing."
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn
        (
            "void BlockDiagGaussSeidelPrecon<Type>::preconditionT\n"
            "(\n"
            "    Field<Type>& xT,\n"
            "    const Field<Type>& bT\n"
            ") const"
         )  << "cannot solve incomplete matrix, no diagonal"
            << abort(FatalError);
    }
}


// ************************************************************************* //
