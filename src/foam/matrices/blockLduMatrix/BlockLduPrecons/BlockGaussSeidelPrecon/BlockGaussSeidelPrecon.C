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

Description
    Gauss-Seidel sweep as a preconditioner.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "BlockGaussSeidelPrecon.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::BlockGaussSeidelPrecon<Type>::calcInvDiag()
{
    typedef CoeffField<Type> TypeCoeffField;
    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::linearType linearType;

    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    // Note: Cannot use inv since the square coefficient requires
    // special treatment.  HJ, 20/Aug/2015

    // Get reference to diagonal
    const TypeCoeffField& d = this->matrix_.diag();

    if (d.activeType() == blockCoeffBase::SCALAR)
    {
        invDiag_.asScalar() = 1/d.asScalar();
    }
    else if (d.activeType() == blockCoeffBase::LINEAR)
    {
        invDiag_.asLinear() =
            cmptDivide
            (
                linearTypeField(d.size(), pTraits<linearType>::one),
                d.asLinear()
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

        linearTypeField lf(activeDiag.size());

        // Take out the diagonal from the diag as a linear type
        contractLinear(lf, activeDiag);

        // Expand diagonal only to full square type and store into luDiag
        expandLinear(luDiag, lf);

        // Keep only off-diagonal in luDiag.
        // Note change of sign to avoid multiplication with -1 when moving
        // to the other side.  HJ, 20/Aug/2015
        luDiag -= activeDiag;

        // Inverse is the inverse of diagonal only
        invDiag_.asLinear() =
            cmptDivide
            (
                linearTypeField(lf.size(), pTraits<linearType>::one),
                lf
            );
    }
    else
    {
        FatalErrorIn
        (
            "void BlockGaussSeidelPrecon<Type>::calcInvDiag()"
        )   << "Problem with coefficient type morphing."
            << abort(FatalError);
    }
}


// Block sweep, symmetric matrix
template<class Type>
template<class DiagType, class ULType>
void Foam::BlockGaussSeidelPrecon<Type>::BlockSweep
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
    // iteration style can be seen in GaussSeidelSolver.C

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
void Foam::BlockGaussSeidelPrecon<Type>::BlockSweep
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
    // iteration style can be seen in GaussSeidelSolver.C

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
Foam::BlockGaussSeidelPrecon<Type>::BlockGaussSeidelPrecon
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
Foam::BlockGaussSeidelPrecon<Type>::BlockGaussSeidelPrecon
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
void Foam::BlockGaussSeidelPrecon<Type>::precondition
(
    Field<Type>& x,
    const Field<Type>& b
) const
{
    typedef CoeffField<Type> TypeCoeffField;

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
                    BlockSweep
                    (
                        x,
                        invDiag_.asScalar(),
                        UpperCoeff.asScalar(),
                        b
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    BlockSweep
                    (
                        x,
                        invDiag_.asScalar(),
                        UpperCoeff.asLinear(),
                        b
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    BlockSweep
                    (
                        x,
                        invDiag_.asScalar(),
                        UpperCoeff.asSquare(),
                        b
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
                    BlockSweep
                    (
                        x,
                        invDiag_.asLinear(),
                        UpperCoeff.asScalar(),
                        b
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    BlockSweep
                    (
                        x,
                        invDiag_.asLinear(),
                        UpperCoeff.asLinear(),
                        b
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    BlockSweep
                    (
                        x,
                        invDiag_.asLinear(),
                        UpperCoeff.asSquare(),
                        b
                    );
                }
            }
        }
        else if (DiagCoeff.activeType() == blockCoeffBase::SQUARE)
        {
            // Add diag coupling to b
            if (bPlusLU_.empty())
            {
                bPlusLU_.setSize(b.size(), pTraits<Type>::zero);
            }

            // Multiply overwrites bPlusLU_: no need to initialise
            // Change of sign accounted via change of sign of bPlusLU_
            // HJ, 20/Aug/2015
            multiply(bPlusLU_, LUDiag_, x);
            bPlusLU_ += b;

            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
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
                    BlockSweep
                    (
                        x,
                        invDiag_.asScalar(),
                        LowerCoeff.asScalar(),
                        UpperCoeff.asScalar(),
                        b
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    BlockSweep
                    (
                        x,
                        invDiag_.asScalar(),
                        LowerCoeff.asLinear(),
                        UpperCoeff.asLinear(),
                        b
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    BlockSweep
                    (
                        x,
                        invDiag_.asScalar(),
                        LowerCoeff.asSquare(),
                        UpperCoeff.asSquare(),
                        b
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
                    BlockSweep
                    (
                        x,
                        invDiag_.asLinear(),
                        LowerCoeff.asScalar(),
                        UpperCoeff.asScalar(),
                        b
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    BlockSweep
                    (
                        x,
                        invDiag_.asLinear(),
                        LowerCoeff.asLinear(),
                        UpperCoeff.asLinear(),
                        b
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                BlockSweep
                (
                    x,
                    invDiag_.asLinear(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    b
                );
            }
        }
        else if (DiagCoeff.activeType() == blockCoeffBase::SQUARE)
        {
            // Add diag coupling to b
            if (bPlusLU_.empty())
            {
                bPlusLU_.setSize(b.size(), pTraits<Type>::zero);
            }

            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
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
//                 Info<< "Diag: " << DiagCoeff.asSquare()[50] << nl
//                     << "LUDiag: " << LUDiag_.asSquare()[50] << nl
//                     << "invDiag: " << invDiag_.asLinear()[50] << nl
//                     << endl;
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
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
                "void BlockGaussSeidelPrecon<Type>::precondition\n"
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
            "void BlockGaussSeidelPrecon<Type>::precondition\n"
            "(\n"
            "    Field<Type>& x,\n"
            "    const Field<Type>& b\n"
            ") const"
         )  << "cannot solve incomplete matrix, no diagonal"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::BlockGaussSeidelPrecon<Type>::preconditionT
(
    Field<Type>& xT,
    const Field<Type>& bT
) const
{
    typedef CoeffField<Type> TypeCoeffField;

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
                    // Transpose multiplication - swap lower and upper coeff
                    BlockSweep
                    (
                        xT,
                        invDiag_.asScalar(),
                        UpperCoeff.asScalar(),
                        LowerCoeff.asScalar(),
                        bT
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Transpose multiplication - swap lower and upper coeff
                    BlockSweep
                    (
                        xT,
                        invDiag_.asScalar(),
                        UpperCoeff.asLinear(),
                        LowerCoeff.asLinear(),
                        bT
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Transpose multiplication - swap lower and upper coeff
                    BlockSweep
                    (
                        xT,
                        invDiag_.asScalar(),
                        UpperCoeff.asSquare(),
                        LowerCoeff.asSquare(),
                        bT
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
                    // Transpose multiplication - swap lower and upper coeff
                    BlockSweep
                    (
                        xT,
                        invDiag_.asLinear(),
                        UpperCoeff.asScalar(),
                        LowerCoeff.asScalar(),
                        bT
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Transpose multiplication - swap lower and upper coeff
                    BlockSweep
                    (
                        xT,
                        invDiag_.asLinear(),
                        UpperCoeff.asLinear(),
                        LowerCoeff.asLinear(),
                        bT
                    );
                }
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    // Transpose multiplication - swap lower and upper coeff
                    BlockSweep
                    (
                        xT,
                        invDiag_.asLinear(),
                        UpperCoeff.asSquare(),
                        LowerCoeff.asSquare(),
                        bT
                    );
                }
            }
        }
        else if (DiagCoeff.activeType() == blockCoeffBase::SQUARE)
        {
            // Add diag coupling to b
            if (bPlusLU_.empty())
            {
                bPlusLU_.setSize(bT.size(), pTraits<Type>::zero);
            }

            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
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
                "void BlockGaussSeidelPrecon<Type>::preconditionT\n"
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
            "void BlockGaussSeidelPrecon<Type>::preconditionT\n"
            "(\n"
            "    Field<Type>& xT,\n"
            "    const Field<Type>& bT\n"
            ") const"
         )  << "cannot solve incomplete matrix, no diagonal"
            << abort(FatalError);
    }
}


// ************************************************************************* //
