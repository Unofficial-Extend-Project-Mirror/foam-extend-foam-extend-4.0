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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "BlockDiagCholeskyPrecon.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::BlockDiagCholeskyPrecon<Type>::calcPreconDiag()
{
    // Note: Assuming lower and upper triangle have the same active type

    typedef CoeffField<Type> TypeCoeffField;

    if (this->matrix_.symmetric())
    {
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                diagMultiply
                (
                    preconDiag_.asScalar(),
                    UpperCoeff.asScalar()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                diagMultiply
                (
                    preconDiag_.asLinear(),
                    UpperCoeff.asLinear()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Transpose multiplication
                diagMultiplyCoeffT
                (
                    preconDiag_.asSquare(),
                    UpperCoeff.asSquare()
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                diagMultiply
                (
                    preconDiag_.asLinear(),
                    UpperCoeff.asScalar()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                diagMultiply
                (
                    preconDiag_.asLinear(),
                    UpperCoeff.asLinear()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Transpose multiplication
                diagMultiplyCoeffT
                (
                    preconDiag_.asSquare(),
                    UpperCoeff.asSquare()
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::SQUARE)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    UpperCoeff.asScalar()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    UpperCoeff.asLinear()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Transpose multiplication
                diagMultiplyCoeffT
                (
                    preconDiag_.asSquare(),
                    UpperCoeff.asSquare()
                );
            }
        }
    }
    else // Asymmetric matrix
    {
        const TypeCoeffField& LowerCoeff = this->matrix_.lower();
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                diagMultiply
                (
                    preconDiag_.asScalar(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                diagMultiply
                (
                    preconDiag_.asLinear(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare()
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                diagMultiply
                (
                    preconDiag_.asLinear(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                diagMultiply
                (
                    preconDiag_.asLinear(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare()
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::SQUARE)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare()
                );
            }
        }
    }

    // Invert the diagonal
//     if (preconDiag_.activeType() == blockCoeffBase::SQUARE)
//     {
//         typedef typename TypeCoeffField::linearTypeField linearTypeField;
//         typedef typename TypeCoeffField::linearType valueType;

//         typedef typename TypeCoeffField::squareTypeField squareTypeField;

//         // Special practice: invert only diagonal of diag coefficient
//         Info<< "Special preconDiag" << endl;

//         // Get reference to active diag
//         const squareTypeField& activeDiag = preconDiag_.asSquare();

//         // Get reference to LU: remove diagonal from active diag
//         squareTypeField& luDiag = LUDiag_.asSquare();

//         linearTypeField lf(preconDiag_.size());

//         // Take out the diagonal from the diag as a linear type
//         contractLinear(lf, activeDiag);

//         expandLinear(luDiag, lf);

//         // Keep only off-diagonal in luDiag
//         // Note change of sign to avoid multiplication with -1 when moving
//         // to the other side.  HJ, 20/Aug/2015
//         luDiag -= activeDiag;

//         // Store inverse of diagonal
//         preconDiag_.clear();

//         // Invert the diagonal part into lf
//         preconDiag_.asLinear() =
//             cmptDivide
//             (
//                 linearTypeField(lf.size(), pTraits<valueType>::one),
//                 lf
//             );
//     }
//     else
    {
        preconDiag_ = inv(preconDiag_);
        LUDiag_.asScalar() = 0;
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::BlockDiagCholeskyPrecon<Type>::diagMultiply
(
    Field<DiagType>& dDiag,
    const Field<ULType>& upper
)
{
    // Precondition the diagonal

    // Get addressing
    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll (upper, coeffI)
    {
        dDiag[upperAddr[coeffI]] -=
            mult.tripleProduct
            (
                upper[coeffI],
                dDiag[lowerAddr[coeffI]],
                upper[coeffI]
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::BlockDiagCholeskyPrecon<Type>::diagMultiplyCoeffT
(
    Field<DiagType>& dDiag,
    const Field<ULType>& upper
)
{
    // Precondition the diagonal

    // Get addressing
    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll (upper, coeffI)
    {
        dDiag[upperAddr[coeffI]] -=
            mult.tripleProduct
            (
                upper[coeffI].T(),        // Upper coefficient transposed
                dDiag[lowerAddr[coeffI]],
                upper[coeffI]
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::BlockDiagCholeskyPrecon<Type>::diagMultiply
(
    Field<DiagType>& dDiag,
    const Field<ULType>& lower,
    const Field<ULType>& upper
)
{
    // Precondition the diagonal

    // Get addressing
    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll (upper, coeffI)
    {
        dDiag[upperAddr[coeffI]] -=
            mult.tripleProduct
            (
                lower[coeffI],
                dDiag[lowerAddr[coeffI]],
                upper[coeffI]
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::BlockDiagCholeskyPrecon<Type>::ILUmultiply
(
    Field<Type>& x,
    const Field<DiagType>& dDiag,
    const Field<ULType>& upper,
    const Field<Type>& b
) const
{
    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll(x, i)
    {
        x[i] = mult(dDiag[i], b[i]);
    }

    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();

    forAll (upper, coeffI)
    {
        x[upperAddr[coeffI]] -=
            mult
            (
                dDiag[upperAddr[coeffI]],
                mult(upper[coeffI], x[lowerAddr[coeffI]])
            );
    }

    forAllReverse (upper, coeffI)
    {
        x[lowerAddr[coeffI]] -=
            mult
            (
                dDiag[lowerAddr[coeffI]],
                mult(upper[coeffI], x[upperAddr[coeffI]])
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::BlockDiagCholeskyPrecon<Type>::ILUmultiplyCoeffT
(
    Field<Type>& x,
    const Field<DiagType>& dDiag,
    const Field<ULType>& upper,
    const Field<Type>& b
) const
{
    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll(x, i)
    {
        x[i] = mult(dDiag[i], b[i]);
    }

    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();

    forAll (upper, coeffI)
    {
        x[upperAddr[coeffI]] -=
            mult
            (
                dDiag[upperAddr[coeffI]],
                mult(upper[coeffI].T(), x[lowerAddr[coeffI]])
            );
    }

    forAllReverse (upper, coeffI)
    {
        x[lowerAddr[coeffI]] -=
            mult
            (
                dDiag[lowerAddr[coeffI]],
                mult(upper[coeffI], x[upperAddr[coeffI]])
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::BlockDiagCholeskyPrecon<Type>::ILUmultiply
(
    Field<Type>& x,
    const Field<DiagType>& dDiag,
    const Field<ULType>& lower,
    const Field<ULType>& upper,
    const Field<Type>& b
) const
{
    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll(x, i)
    {
        x[i] = mult(dDiag[i], b[i]);
    }

    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();
    const unallocLabelList& losortAddr = this->matrix_.lduAddr().losortAddr();

    register label losortCoeff;

    forAll (lower, coeffI)
    {
        losortCoeff = losortAddr[coeffI];

        x[upperAddr[losortCoeff]] -=
            mult
            (
                dDiag[upperAddr[losortCoeff]],
                mult(lower[losortCoeff], x[lowerAddr[losortCoeff]])
            );
    }

    forAllReverse (upper, coeffI)
    {
        x[lowerAddr[coeffI]] -=
            mult
            (
                dDiag[lowerAddr[coeffI]],
                mult(upper[coeffI], x[upperAddr[coeffI]])
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::BlockDiagCholeskyPrecon<Type>::ILUmultiplyTranspose
(
    Field<Type>& x,
    const Field<DiagType>& dDiag,
    const Field<ULType>& lower,
    const Field<ULType>& upper,
    const Field<Type>& b
) const
{
    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll(x, i)
    {
        x[i] = mult(dDiag[i], b[i]);
    }

    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();
    const unallocLabelList& losortAddr = this->matrix_.lduAddr().losortAddr();

    register label losortCoeff;

    //HJ Not sure if the coefficient itself needs to be transposed.
    // HJ, 30/Oct/2007
    forAll (lower, coeffI)
    {
        x[upperAddr[coeffI]] -=
            mult
            (
                dDiag[upperAddr[coeffI]],
                mult(upper[coeffI], x[lowerAddr[coeffI]])
            );
    }

    forAllReverse (upper, coeffI)
    {
        losortCoeff = losortAddr[coeffI];

        x[lowerAddr[losortCoeff]] -=
            mult
            (
                dDiag[lowerAddr[losortCoeff]],
                mult(lower[losortCoeff], x[upperAddr[losortCoeff]])
            );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockDiagCholeskyPrecon<Type>::BlockDiagCholeskyPrecon
(
    const BlockLduMatrix<Type>& matrix
)
:
    BlockLduPrecon<Type>(matrix),
    preconDiag_(matrix.diag()),
    LUDiag_(matrix.lduAddr().size()),
    bPlusLU_()
{
    calcPreconDiag();
}


template<class Type>
Foam::BlockDiagCholeskyPrecon<Type>::BlockDiagCholeskyPrecon
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
:
    BlockLduPrecon<Type>(matrix),
    preconDiag_(matrix.diag()),
    LUDiag_(matrix.lduAddr().size()),
    bPlusLU_()
{
    calcPreconDiag();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockDiagCholeskyPrecon<Type>::~BlockDiagCholeskyPrecon()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::BlockDiagCholeskyPrecon<Type>::precondition
(
    Field<Type>& x,
    const Field<Type>& b
) const
{
    typedef CoeffField<Type> TypeCoeffField;

    // Note: Assuming lower and upper triangle have the same active type

    if (this->matrix_.symmetric())
    {
        const TypeCoeffField& DiagCoeff = this->matrix_.diag();
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        if (DiagCoeff.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asScalar(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asScalar(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Multiplication using transposed upper square coefficient
                ILUmultiplyCoeffT
                (
                    x,
                    preconDiag_.asScalar(),
                    UpperCoeff.asSquare(),
                    b
                );
            }
        }
        else if (DiagCoeff.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Multiplication using transposed upper square coefficient
                ILUmultiplyCoeffT
                (
                    x,
                    preconDiag_.asLinear(),
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

            // Multiply overwrites bPlusLU_: no need to initialise
            // Change of sign accounted via change of sign of bPlusLU_
            // HJ, 20/Aug/2015
            multiply(bPlusLU_, LUDiag_, x);
            bPlusLU_ += b;

            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                // Note linear preconDiag due to decoupling
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    UpperCoeff.asScalar(),
                    bPlusLU_
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                // Note linear preconDiag due to decoupling
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    UpperCoeff.asLinear(),
                    bPlusLU_
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Multiplication using transposed upper square coefficient
                // Note linear preconDiag due to decoupling
                ILUmultiplyCoeffT
                (
                    x,
                    preconDiag_.asLinear(),
                    UpperCoeff.asSquare(),
                    bPlusLU_
                );
            }
        }
    }
    else // Asymmetric matrix
    {
        const TypeCoeffField& DiagCoeff = this->matrix_.diag();
        const TypeCoeffField& LowerCoeff = this->matrix_.lower();
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        if (DiagCoeff.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asScalar(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asScalar(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asScalar(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    b
                );
            }
        }
        else if (DiagCoeff.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
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

            // Multiply overwrites bPlusLU_: no need to initialise
            // Change of sign accounted via change of sign of bPlusLU_
            // HJ, 20/Aug/2015
            multiply(bPlusLU_, LUDiag_, x);
            bPlusLU_ += b;

            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                // Note linear preconDiag due to decoupling
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    bPlusLU_
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                // Note linear preconDiag due to decoupling
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    bPlusLU_
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Note linear preconDiag due to decoupling
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    bPlusLU_
                );
            }
        }
    }
}


template<class Type>
void Foam::BlockDiagCholeskyPrecon<Type>::preconditionT
(
    Field<Type>& xT,
    const Field<Type>& bT
) const
{
    typedef CoeffField<Type> TypeCoeffField;

    // Note: Assuming lower and upper triangle have the same active type

    if (this->matrix_.symmetric())
    {
        precondition(xT, bT);
    }
    else // Asymmetric matrix
    {
        const TypeCoeffField& DiagCoeff = this->matrix_.diag();
        const TypeCoeffField& LowerCoeff = this->matrix_.lower();
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        if (DiagCoeff.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asScalar(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asScalar(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asScalar(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    bT
                );
            }
        }
        else if (DiagCoeff.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asLinear(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asLinear(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asLinear(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    bT
                );
            }
        }
        else if (DiagCoeff.activeType() == blockCoeffBase::SQUARE)
        {
            // Add diag coupling to b
            if (bPlusLU_.empty())
            {
                bPlusLU_.setSize(bT.size(), pTraits<Type>::zero);
            }

            // Multiply overwrites bPlusLU_: no need to initialise
            // Change of sign accounted via change of sign of bPlusLU_
            // HJ, 20/Aug/2015
            multiply(bPlusLU_, LUDiag_, xT);
            bPlusLU_ += bT;

            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                // Note linear preconDiag due to decoupling
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asLinear(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    bPlusLU_
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                // Note linear preconDiag due to decoupling
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asLinear(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    bPlusLU_
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Note linear preconDiag due to decoupling
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asLinear(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    bPlusLU_
                );
            }
        }
    }
}


// ************************************************************************* //
