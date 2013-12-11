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
    Gauss-Seidel sweep as a preconditioner.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "BlockGaussSeidelPrecon.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

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
    const unallocLabelList& ownStart = this->matrix_.lduAddr().ownerStartAddr();

    const label nRows = ownStart.size() - 1;

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    // Klas Jareteg: 2013-02-10:
    // Must transfer data between the different CPUs. Notes on the Jacobi
    // iteration style can be seen in GaussSeidelSolver.C

    for (label sweep = 0; sweep < nSweeps_; sweep++)
    {
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

            // Distribute the neighbour side using current x
            for (curCoeff = fStart; curCoeff < fEnd; curCoeff++)
            {
                // lower = upper transposed
                bPrime_[u[curCoeff]] -=
                    mult(mult.transpose(upper[curCoeff]), curX);
            }
        }
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
    const unallocLabelList& ownStart = this->matrix_.lduAddr().ownerStartAddr();

    const label nRows = ownStart.size() - 1;

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    // Klas Jareteg: 2013-02-10:
    // Must transfer data between the different CPUs. Notes on the Jacobi
    // iteration style can be seen in GaussSeidelSolver.C

    for (label sweep = 0; sweep < nSweeps_; sweep++)
    {
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

            // Distribute the neighbour side using current x
            for (curCoeff = fStart; curCoeff < fEnd; curCoeff++)
            {
                bPrime_[u[curCoeff]] -= mult(lower[curCoeff], curX);
            }
        }
    }
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
        TypeCoeffField dDCoeff = inv(this->matrix_.diag());

        multiply(x, dDCoeff, b);
    }
    else if (this->matrix_.symmetric())
    {
        TypeCoeffField dDCoeff = inv(this->matrix_.diag());

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

        if (dDCoeff.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                BlockSweep
                (
                    x,
                    dDCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                BlockSweep
                (
                    x,
                    dDCoeff.asScalar(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                BlockSweep
                (
                    x,
                    dDCoeff.asScalar(),
                    UpperCoeff.asSquare(),
                    b
                );
            }
        }
        else if (dDCoeff.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                BlockSweep
                (
                    x,
                    dDCoeff.asLinear(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                BlockSweep
                (
                    x,
                    dDCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                BlockSweep
                (
                    x,
                    dDCoeff.asLinear(),
                    UpperCoeff.asSquare(),
                    b
                );
            }
        }
        else if (dDCoeff.activeType() == blockCoeffBase::SQUARE)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                BlockSweep
                (
                    x,
                    dDCoeff.asSquare(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                BlockSweep
                (
                    x,
                    dDCoeff.asSquare(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                BlockSweep
                (
                    x,
                    dDCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    b
                );
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
        TypeCoeffField dDCoeff = inv(this->matrix_.diag());

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

        if (dDCoeff.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                BlockSweep
                (
                    x,
                    dDCoeff.asScalar(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                BlockSweep
                (
                    x,
                    dDCoeff.asScalar(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                BlockSweep
                (
                    x,
                    dDCoeff.asScalar(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    b
                );
            }
        }
        else if (dDCoeff.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                BlockSweep
                (
                    x,
                    dDCoeff.asLinear(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                BlockSweep
                (
                    x,
                    dDCoeff.asLinear(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                BlockSweep
                (
                    x,
                    dDCoeff.asLinear(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    b
                );
            }
        }
        else if (dDCoeff.activeType() == blockCoeffBase::SQUARE)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                BlockSweep
                (
                    x,
                    dDCoeff.asSquare(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                BlockSweep
                (
                    x,
                    dDCoeff.asSquare(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                BlockSweep
                (
                    x,
                    dDCoeff.asSquare(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    b
                );
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
        TypeCoeffField dDCoeff = inv(this->matrix_.diag());

        multiply(xT, dDCoeff, bT);
    }
    else if (this->matrix_.symmetric() || this->matrix_.asymmetric())
    {
        TypeCoeffField dDCoeff = inv(this->matrix_.diag());

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

        if (dDCoeff.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                // Transpose multiplication - swap lower and upper coeff arrays
                BlockSweep
                (
                    xT,
                    dDCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    LowerCoeff.asScalar(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                // Transpose multiplication - swap lower and upper coeff arrays
                BlockSweep
                (
                    xT,
                    dDCoeff.asScalar(),
                    UpperCoeff.asLinear(),
                    LowerCoeff.asLinear(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Transpose multiplication - swap lower and upper coeff arrays
                BlockSweep
                (
                    xT,
                    dDCoeff.asScalar(),
                    UpperCoeff.asSquare(),
                    LowerCoeff.asSquare(),
                    bT
                );
            }
        }
        else if (dDCoeff.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                // Transpose multiplication - swap lower and upper coeff arrays
                BlockSweep
                (
                    xT,
                    dDCoeff.asLinear(),
                    UpperCoeff.asScalar(),
                    LowerCoeff.asScalar(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                // Transpose multiplication - swap lower and upper coeff arrays
                BlockSweep
                (
                    xT,
                    dDCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    LowerCoeff.asLinear(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Transpose multiplication - swap lower and upper coeff arrays
                BlockSweep
                (
                    xT,
                    dDCoeff.asLinear(),
                    UpperCoeff.asSquare(),
                    LowerCoeff.asSquare(),
                    bT
                );
            }
        }
        else if (dDCoeff.activeType() == blockCoeffBase::SQUARE)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                // Transpose multiplication - swap lower and upper coeff arrays
                BlockSweep
                (
                    xT,
                    dDCoeff.asSquare(),
                    UpperCoeff.asScalar(),
                    LowerCoeff.asScalar(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                // Transpose multiplication - swap lower and upper coeff arrays
                BlockSweep
                (
                    xT,
                    dDCoeff.asSquare(),
                    UpperCoeff.asLinear(),
                    LowerCoeff.asLinear(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Transpose multiplication - swap lower and upper coeff arrays
                BlockSweep
                (
                    xT,
                    dDCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    LowerCoeff.asSquare(),
                    bT
                );
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
