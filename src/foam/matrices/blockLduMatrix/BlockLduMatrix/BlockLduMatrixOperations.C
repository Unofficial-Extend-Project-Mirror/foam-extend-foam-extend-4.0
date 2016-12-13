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

\*---------------------------------------------------------------------------*/

#include "BlockLduMatrix.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::BlockLduMatrix<Type>::sumDiag()
{
    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    TypeCoeffField& Diag = this->diag();

    const unallocLabelList& l = lduAddr().lowerAddr();
    const unallocLabelList& u = lduAddr().upperAddr();

    if (this->symmetric())
    {
        // Symmetric matrix: re-use upper transpose for lower coefficients

        const TypeCoeffField& Upper =
            const_cast<const BlockLduMatrix<Type>&>(*this).upper();

        if
        (
            Upper.activeType() == blockCoeffBase::SQUARE
         || Diag.activeType() == blockCoeffBase::SQUARE
        )
        {
            const squareTypeField& activeUpper = Upper.asSquare();
            squareTypeField& activeDiag = Diag.asSquare();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeDiag[l[coeffI]] += activeUpper[coeffI].T();
                activeDiag[u[coeffI]] += activeUpper[coeffI];
            }
        }
        else if
        (
            Upper.activeType() == blockCoeffBase::LINEAR
         || Diag.activeType() == blockCoeffBase::LINEAR
        )
        {
            const linearTypeField& activeUpper = Upper.asLinear();
            linearTypeField& activeDiag = Diag.asLinear();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeDiag[l[coeffI]] += activeUpper[coeffI];
                activeDiag[u[coeffI]] += activeUpper[coeffI];
            }
        }
        else if
        (
            Upper.activeType() == blockCoeffBase::SCALAR
         || Diag.activeType() == blockCoeffBase::SCALAR
        )
        {
            const scalarTypeField& activeUpper = Upper.asScalar();
            scalarTypeField& activeDiag = Diag.asScalar();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeDiag[l[coeffI]] += activeUpper[coeffI];
                activeDiag[u[coeffI]] += activeUpper[coeffI];
            }
        }
    }
    else if (this->asymmetric())
    {
        // Full asymmetric matrix

        const TypeCoeffField& Lower =
            const_cast<const BlockLduMatrix<Type>&>(*this).lower();

        const TypeCoeffField& Upper =
            const_cast<const BlockLduMatrix<Type>&>(*this).upper();

        if
        (
            Lower.activeType() == blockCoeffBase::SQUARE
         || Upper.activeType() == blockCoeffBase::SQUARE
         || Diag.activeType() == blockCoeffBase::SQUARE
        )
        {
            const squareTypeField& activeLower = Lower.asSquare();
            const squareTypeField& activeUpper = Upper.asSquare();
            squareTypeField& activeDiag = Diag.asSquare();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeDiag[l[coeffI]] += activeLower[coeffI];
                activeDiag[u[coeffI]] += activeUpper[coeffI];
            }
        }
        else if
        (
            Lower.activeType() == blockCoeffBase::LINEAR
         || Upper.activeType() == blockCoeffBase::LINEAR
         || Diag.activeType() == blockCoeffBase::LINEAR
        )
        {
            const linearTypeField& activeLower = Lower.asLinear();
            const linearTypeField& activeUpper = Upper.asLinear();
            linearTypeField& activeDiag = Diag.asLinear();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeDiag[l[coeffI]] += activeLower[coeffI];
                activeDiag[u[coeffI]] += activeUpper[coeffI];
            }
        }
        else if
        (
            Lower.activeType() == blockCoeffBase::SCALAR
         || Upper.activeType() == blockCoeffBase::SCALAR
         || Diag.activeType() == blockCoeffBase::SCALAR
        )
        {
            const scalarTypeField& activeLower = Lower.asScalar();
            const scalarTypeField& activeUpper = Upper.asScalar();
            scalarTypeField& activeDiag = Diag.asScalar();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeDiag[l[coeffI]] += activeLower[coeffI];
                activeDiag[u[coeffI]] += activeUpper[coeffI];
            }
        }
    }
    else
    {
        FatalErrorIn("void BlockLduMatrix<Type>::sumDiag()")
            << "No off-diagonal available"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::BlockLduMatrix<Type>::negSumDiag()
{
    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    TypeCoeffField& Diag = this->diag();

    const unallocLabelList& l = lduAddr().lowerAddr();
    const unallocLabelList& u = lduAddr().upperAddr();

    if (this->symmetric())
    {
        // Symmetric matrix: re-use upper transpose for lower coefficients

        const TypeCoeffField& Upper =
            const_cast<const BlockLduMatrix<Type>&>(*this).upper();

        if
        (
            Upper.activeType() == blockCoeffBase::SQUARE
         || Diag.activeType() == blockCoeffBase::SQUARE
        )
        {
            const squareTypeField& activeUpper = Upper.asSquare();
            squareTypeField& activeDiag = Diag.asSquare();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeDiag[l[coeffI]] -= activeUpper[coeffI].T();
                activeDiag[u[coeffI]] -= activeUpper[coeffI];
            }
        }
        else if
        (
            Upper.activeType() == blockCoeffBase::LINEAR
         || Diag.activeType() == blockCoeffBase::LINEAR
        )
        {
            const linearTypeField& activeUpper = Upper.asLinear();
            linearTypeField& activeDiag = Diag.asLinear();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeDiag[l[coeffI]] -= activeUpper[coeffI];
                activeDiag[u[coeffI]] -= activeUpper[coeffI];
            }
        }
        else if
        (
            Upper.activeType() == blockCoeffBase::SCALAR
         || Diag.activeType() == blockCoeffBase::SCALAR
        )
        {
            const scalarTypeField& activeUpper = Upper.asScalar();
            scalarTypeField& activeDiag = Diag.asScalar();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeDiag[l[coeffI]] -= activeUpper[coeffI];
                activeDiag[u[coeffI]] -= activeUpper[coeffI];
            }
        }
    }
    else if (this->asymmetric())
    {
        // Full asymmetric matrix

        const TypeCoeffField& Lower =
            const_cast<const BlockLduMatrix<Type>&>(*this).lower();

        const TypeCoeffField& Upper =
            const_cast<const BlockLduMatrix<Type>&>(*this).upper();

        if
        (
            Lower.activeType() == blockCoeffBase::SQUARE
         || Upper.activeType() == blockCoeffBase::SQUARE
         || Diag.activeType() == blockCoeffBase::SQUARE
        )
        {
            const squareTypeField& activeLower = Lower.asSquare();
            const squareTypeField& activeUpper = Upper.asSquare();
            squareTypeField& activeDiag = Diag.asSquare();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeDiag[l[coeffI]] -= activeLower[coeffI];
                activeDiag[u[coeffI]] -= activeUpper[coeffI];
            }
        }
        else if
        (
            Lower.activeType() == blockCoeffBase::LINEAR
         || Upper.activeType() == blockCoeffBase::LINEAR
         || Diag.activeType() == blockCoeffBase::LINEAR
        )
        {
            const linearTypeField& activeLower = Lower.asLinear();
            const linearTypeField& activeUpper = Upper.asLinear();
            linearTypeField& activeDiag = Diag.asLinear();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeDiag[l[coeffI]] -= activeLower[coeffI];
                activeDiag[u[coeffI]] -= activeUpper[coeffI];
            }
        }
        else if
        (
            Lower.activeType() == blockCoeffBase::SCALAR
         || Upper.activeType() == blockCoeffBase::SCALAR
         || Diag.activeType() == blockCoeffBase::SCALAR
        )
        {
            const scalarTypeField& activeLower = Lower.asScalar();
            const scalarTypeField& activeUpper = Upper.asScalar();
            scalarTypeField& activeDiag = Diag.asScalar();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeDiag[l[coeffI]] -= activeLower[coeffI];
                activeDiag[u[coeffI]] -= activeUpper[coeffI];
            }
        }
    }
    else
    {
        FatalErrorIn("void BlockLduMatrix<Type>::negSumDiag()")
            << "No off-diagonal available"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::BlockLduMatrix<Type>::check() const
{
    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    // TODO: Account for coupled boundary patch coeffs
    // HJ, 29/Oct/2015

    // Get original diag
    const TypeCoeffField& d = this->diag();

    // Copy the diagonal.  It is initialised to zero
    TypeCoeffField SumOffDiag(d.size());

    const unallocLabelList& l = lduAddr().lowerAddr();
    const unallocLabelList& u = lduAddr().upperAddr();

    if (this->symmetric())
    {
        // Symmetric matrix: re-use upper transpose for lower coefficients

        const TypeCoeffField& Upper = this->upper();

        if
        (
            Upper.activeType() == blockCoeffBase::SQUARE
         || d.activeType() == blockCoeffBase::SQUARE
        )
        {
            // Note: For a square coefficient, the matrix needs to be analysed
            // row-by-row, including the contribution of the off-diagonal
            // elements in the diagonal matrix

            // Get result as linear
            linearTypeField& activeSumOffDiag = SumOffDiag.asLinear();

            // Do diagonal elements

            // Create the diagonal matrix element without the diagonal
            const squareTypeField& activeDiag = d.asSquare();

            linearTypeField diagDiag(activeDiag.size());
            squareTypeField luDiag(activeDiag.size());

            // Take out the diagonal from the diag as a linear type
            contractLinear(diagDiag, activeDiag);

            // Expand diagonal only to full square type and store into luDiag
            expandLinear(luDiag, diagDiag);

            // Keep only off-diagonal in luDiag
            luDiag = activeDiag - luDiag;

            sumMagToDiag(activeSumOffDiag, luDiag);

            // Do the off-diagonal elements by collapsing each row
            // into the linear form

            const squareTypeField& activeUpper = Upper.asSquare();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeSumOffDiag[l[coeffI]] +=
                    sumMagToDiag(activeUpper[coeffI].T());

                activeSumOffDiag[u[coeffI]] +=
                    sumMagToDiag(activeUpper[coeffI]);
            }

            // Check diagonal dominance

            diagDiag = cmptMag(diagDiag);

            // Divide diagonal with sum of off-diagonals
            cmptDivide(diagDiag, diagDiag, activeSumOffDiag);

            InfoIn("void BlockLduMatrix<Type>::check() const)")
                << "Symmetric matrix: " << activeDiag.size()
                << " diagonal dominance sym square: "
                << Foam::min(diagDiag)
                << endl;
        }
        else if
        (
            Upper.activeType() == blockCoeffBase::LINEAR
         || d.activeType() == blockCoeffBase::LINEAR
        )
        {
            const linearTypeField& activeDiag = d.asLinear();

            const linearTypeField& activeUpper = Upper.asLinear();
            linearTypeField& activeSumOffDiag = SumOffDiag.asLinear();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeSumOffDiag[l[coeffI]] += cmptMag(activeUpper[coeffI]);
                activeSumOffDiag[u[coeffI]] += cmptMag(activeUpper[coeffI]);
            }

            linearTypeField diagDiag = cmptMag(activeDiag);

            cmptDivide(diagDiag, diagDiag, activeSumOffDiag);

            InfoIn("void BlockLduMatrix<Type>::check() const)")
                << "Symmetric matrix: " << activeDiag.size()
                << " diagonal dominance sym linear: "
                << Foam::min(diagDiag)
                << endl;
        }
        else if
        (
            Upper.activeType() == blockCoeffBase::SCALAR
         || d.activeType() == blockCoeffBase::SCALAR
        )
        {
            const scalarTypeField& activeDiag = d.asScalar();

            const scalarTypeField& activeUpper = Upper.asScalar();
            scalarTypeField& activeSumOffDiag = SumOffDiag.asScalar();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeSumOffDiag[l[coeffI]] += cmptMag(activeUpper[coeffI]);
                activeSumOffDiag[u[coeffI]] += cmptMag(activeUpper[coeffI]);
            }

            InfoIn("void BlockLduMatrix<Type>::check() const)")
                << "Symmetric matrix: " << activeDiag.size()
                << " diagonal dominance sym scalar: "
                << Foam::min(mag(activeDiag)/activeSumOffDiag)
                << endl;
        }
    }
    else if (this->asymmetric())
    {
        // Full asymmetric matrix

        const TypeCoeffField& Lower = this->lower();
        const TypeCoeffField& Upper = this->upper();

        if
        (
            Lower.activeType() == blockCoeffBase::SQUARE
         || Upper.activeType() == blockCoeffBase::SQUARE
         || d.activeType() == blockCoeffBase::SQUARE
        )
        {
            // Note: For a square coefficient, the matrix needs to be analysed
            // row-by-row, including the contribution of the off-diagonal
            // elements in the diagonal matrix

            // Get result as linear
            linearTypeField& activeSumOffDiag = SumOffDiag.asLinear();

            // Do diagonal elements

            // Create the diagonal matrix element without the diagonal
            const squareTypeField& activeDiag = d.asSquare();

            linearTypeField diagDiag(activeDiag.size());
            squareTypeField luDiag(activeDiag.size());

            // Take out the diagonal from the diag as a linear type
            contractLinear(diagDiag, activeDiag);

            // Expand diagonal only to full square type and store into luDiag
            expandLinear(luDiag, diagDiag);

            // Keep only off-diagonal in luDiag
            luDiag = activeDiag - luDiag;

            sumMagToDiag(activeSumOffDiag, luDiag);

            // Do the off-diagonal elements by collapsing each row
            // into the linear form

            const squareTypeField& activeLower = Lower.asSquare();
            const squareTypeField& activeUpper = Upper.asSquare();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeSumOffDiag[l[coeffI]] +=
                    sumMagToDiag(activeLower[coeffI]);

                activeSumOffDiag[u[coeffI]] +=
                    sumMagToDiag(activeUpper[coeffI]);
            }

            // Check diagonal dominance

            diagDiag = cmptMag(diagDiag);

            // Divide diagonal with sum of off-diagonals
            cmptDivide(diagDiag, diagDiag, activeSumOffDiag);

            InfoIn("void BlockLduMatrix<Type>::check() const)")
                << "Asymmetric matrix: " << activeDiag.size()
                << " diagonal dominance assym square: "
                << Foam::min(diagDiag)
                << endl;
        }
        else if
        (
            Lower.activeType() == blockCoeffBase::LINEAR
         || Upper.activeType() == blockCoeffBase::LINEAR
         || d.activeType() == blockCoeffBase::LINEAR
        )
        {
            const linearTypeField& activeDiag = d.asLinear();

            const linearTypeField& activeLower = Lower.asLinear();
            const linearTypeField& activeUpper = Upper.asLinear();
            linearTypeField& activeSumOffDiag = SumOffDiag.asLinear();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeSumOffDiag[l[coeffI]] += cmptMag(activeLower[coeffI]);
                activeSumOffDiag[u[coeffI]] += cmptMag(activeUpper[coeffI]);
            }

            linearTypeField diagDiag = cmptMag(activeDiag);

            cmptDivide(diagDiag, diagDiag, activeSumOffDiag);

            InfoIn("void BlockLduMatrix<Type>::check() const)")
                << "Asymmetric matrix: " << activeDiag.size()
                << " diagonal dominance assym linear: "
                << Foam::min(diagDiag)
                << endl;
        }
        else if
        (
            Lower.activeType() == blockCoeffBase::SCALAR
         || Upper.activeType() == blockCoeffBase::SCALAR
         || d.activeType() == blockCoeffBase::SCALAR
        )
        {
            const scalarTypeField& activeDiag = d.asScalar();

            const scalarTypeField& activeLower = Lower.asScalar();
            const scalarTypeField& activeUpper = Upper.asScalar();
            scalarTypeField& activeSumOffDiag = SumOffDiag.asScalar();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeSumOffDiag[l[coeffI]] += cmptMag(activeLower[coeffI]);
                activeSumOffDiag[u[coeffI]] += cmptMag(activeUpper[coeffI]);
            }

            InfoIn("void BlockLduMatrix<Type>::check() const)")
                << "Asymmetric matrix: "  << activeDiag.size()
                << " diagonal dominance assym scalar: "
                << Foam::min(mag(activeDiag)/activeSumOffDiag)
                << endl;
        }
    }
    else
    {
        InfoIn("void BlockLduMatrix<Type>::check() const)")
            << "Diagonal matrix" << endl;
    }
}


template<class Type>
void Foam::BlockLduMatrix<Type>::relax
(
    const Field<Type>& x,
    Field<Type>& b,
    const scalar alpha
)
{
    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    //HJ Missing code: add coupling coefficients to under-relaxation
    //   HJ, 21/Feb/2008

    if (alpha <= 0)
    {
        return;
    }

    TypeCoeffField& Diag = this->diag();

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    const unallocLabelList& l = lduAddr().lowerAddr();
    const unallocLabelList& u = lduAddr().upperAddr();

    if (this->symmetric())
    {
        // Symmetric matrix: re-use upper transpose for lower coefficients

        const TypeCoeffField& Upper =
            const_cast<const BlockLduMatrix<Type>&>(*this).upper();

        if
        (
            Upper.activeType() == blockCoeffBase::SQUARE
         || Diag.activeType() == blockCoeffBase::SQUARE
        )
        {
            const squareTypeField& activeUpper = Upper.asSquare();
            squareTypeField& activeDiag = Diag.asSquare();

            // Make a copy of diagonal before relaxation
            squareTypeField activeDiagOld = activeDiag;

            // There are options for under-relaxing the block diagonal
            // coefficient.  Currently, the complete diagonal block is
            // under-relaxed.  There's no checking on the off-diag sum

            squareTypeField sumOff
            (
                activeDiag.size(),
                pTraits<typename TypeCoeffField::squareType>::zero
            );

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                sumOff[u[coeffI]] += cmptMag(activeUpper[coeffI].T());
                sumOff[l[coeffI]] += cmptMag(activeUpper[coeffI]);
            }

            // Reconsider under-relaxation of square blocks.
            // HJ, 23/Sep/2011 (2 places)
            activeDiag = max(activeDiag, sumOff);
            activeDiag *= 1.0/alpha;

            // Add the relaxation contribution to b
            forAll (b, i)
            {
                b[i] += mult(activeDiag[i] - activeDiagOld[i], x[i]);
            }
        }
        else if
        (
            Upper.activeType() == blockCoeffBase::LINEAR
         || Diag.activeType() == blockCoeffBase::LINEAR
        )
        {
            const linearTypeField& activeUpper = Upper.asLinear();
            linearTypeField& activeDiag = Diag.asLinear();

            // Make a copy of diagonal before relaxation
            linearTypeField activeDiagOld = activeDiag;

            linearTypeField sumOff
            (
                activeDiag.size(),
                pTraits<typename TypeCoeffField::linearType>::zero
            );

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                sumOff[u[coeffI]] += cmptMag(activeUpper[coeffI]);
                sumOff[l[coeffI]] += cmptMag(activeUpper[coeffI]);
            }

            activeDiag = max(activeDiag, sumOff);
            activeDiag *= 1.0/alpha;

            // Add the relaxation contribution to b
            forAll (b, i)
            {
                b[i] += mult(activeDiag[i] - activeDiagOld[i], x[i]);
            }
        }
        else if
        (
            Upper.activeType() == blockCoeffBase::SCALAR
         || Diag.activeType() == blockCoeffBase::SCALAR
        )
        {
            const scalarTypeField& activeUpper = Upper.asScalar();
            scalarTypeField& activeDiag = Diag.asScalar();

            // Make a copy of diagonal before relaxation
            scalarTypeField activeDiagOld = activeDiag;

            scalarTypeField sumOff
            (
                activeDiag.size(),
                pTraits<typename TypeCoeffField::scalarType>::zero
            );

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                sumOff[u[coeffI]] += mag(activeUpper[coeffI]);
                sumOff[l[coeffI]] += mag(activeUpper[coeffI]);
            }

            activeDiag = max(activeDiag, sumOff);
            activeDiag *= 1.0/alpha;

            // Add the relaxation contribution to b
            forAll (b, i)
            {
                b[i] += (activeDiag[i] - activeDiagOld[i])*x[i];
            }
        }
    }
    else if (this->asymmetric())
    {
        // Full asymmetric matrix

        const TypeCoeffField& Lower =
            const_cast<const BlockLduMatrix<Type>&>(*this).lower();

        const TypeCoeffField& Upper =
            const_cast<const BlockLduMatrix<Type>&>(*this).upper();

        if
        (
            Lower.activeType() == blockCoeffBase::SQUARE
         || Upper.activeType() == blockCoeffBase::SQUARE
         || Diag.activeType() == blockCoeffBase::SQUARE
        )
        {
            const squareTypeField& activeLower = Lower.asSquare();
            const squareTypeField& activeUpper = Upper.asSquare();
            squareTypeField& activeDiag = Diag.asSquare();

            // Make a copy of diagonal before relaxation
            squareTypeField activeDiagOld = activeDiag;

            // There are options for under-relaxing the block diagonal
            // coefficient.  Currently, the complete diagonal block is
            // under-relaxed.  There's no checking on the off-diag sum

            squareTypeField sumOff
            (
                activeDiag.size(),
                pTraits<typename TypeCoeffField::squareType>::zero
            );

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                sumOff[u[coeffI]] += cmptMag(activeLower[coeffI]);
                sumOff[l[coeffI]] += cmptMag(activeUpper[coeffI]);
            }

            // Reconsider under-relaxation of square blocks.
            // HJ, 23/Sep/2011 (2 places)
            activeDiag = max(activeDiag, sumOff);
            activeDiag *= 1.0/alpha;

            // Add the relaxation contribution to b
            forAll (b, i)
            {
                b[i] += mult(activeDiag[i] - activeDiagOld[i], x[i]);
            }
        }
        else if
        (
            Lower.activeType() == blockCoeffBase::LINEAR
         || Upper.activeType() == blockCoeffBase::LINEAR
         || Diag.activeType() == blockCoeffBase::LINEAR
        )
        {
            const linearTypeField& activeLower = Lower.asLinear();
            const linearTypeField& activeUpper = Upper.asLinear();
            linearTypeField& activeDiag = Diag.asLinear();

            // Make a copy of diagonal before relaxation
            linearTypeField activeDiagOld = activeDiag;

            linearTypeField sumOff
            (
                activeDiag.size(),
                pTraits<typename TypeCoeffField::linearType>::zero
            );

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                sumOff[u[coeffI]] += cmptMag(activeLower[coeffI]);
                sumOff[l[coeffI]] += cmptMag(activeUpper[coeffI]);
            }

            activeDiag = max(activeDiag, sumOff);
            activeDiag *= 1.0/alpha;

            // Add the relaxation contribution to b
            forAll (b, i)
            {
                b[i] += mult(activeDiag[i] - activeDiagOld[i], x[i]);
            }
        }
        else if
        (
            Lower.activeType() == blockCoeffBase::SCALAR
         || Upper.activeType() == blockCoeffBase::SCALAR
         || Diag.activeType() == blockCoeffBase::SCALAR
        )
        {
            const scalarTypeField& activeLower = Lower.asScalar();
            const scalarTypeField& activeUpper = Upper.asScalar();
            scalarTypeField& activeDiag = Diag.asScalar();

            // Make a copy of diagonal before relaxation
            scalarTypeField activeDiagOld = activeDiag;

            scalarTypeField sumOff
            (
                activeDiag.size(),
                pTraits<typename TypeCoeffField::scalarType>::zero
            );

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                sumOff[u[coeffI]] += mag(activeLower[coeffI]);
                sumOff[l[coeffI]] += mag(activeUpper[coeffI]);
            }

            activeDiag = max(activeDiag, sumOff);
            activeDiag *= 1.0/alpha;

            // Add the relaxation contribution to b
            forAll (b, i)
            {
                b[i] += activeDiag[i] - activeDiagOld[i]*x[i];
            }
        }
    }
}


template<class Type>
void Foam::BlockLduMatrix<Type>::setValue
(
    const label eqnIndex,
    const Type& value
)
{
    BlockConstraint<Type> cp(eqnIndex, value);

    if (!fixedEqns_.found(eqnIndex))
    {
        fixedEqns_.insert(eqnIndex, cp);
    }
    else
    {
        WarningIn
        (
            "void BlockLduMatrix<Type>::setValue(const label eqnIndex, "
            "const Type& value)"
        )   << "Adding constraint on an already constrained point."
            << "  Point: " << eqnIndex << endl;

        fixedEqns_[eqnIndex].combine(cp);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::BlockLduMatrix<Type>::residual
(
    const Field<Type>& x
) const
{
    Field<Type> Ax(x.size());
    Amul(Ax, x);
    return -Ax;
}


template<class Type>
typename Foam::tmp<Foam::Field<Type> > Foam::BlockLduMatrix<Type>::residual
(
    const Field<Type>& x,
    const Field<Type>& b
) const
{
    return b + residual(x);
}


template<class Type>
void Foam::BlockLduMatrix<Type>::negate()
{
    if (lowerPtr_)
    {
        lowerPtr_->negate();
    }

    if (upperPtr_)
    {
        upperPtr_->negate();
    }

    if (diagPtr_)
    {
        diagPtr_->negate();
    }

    // Do coupling coefficients
    coupleUpper_.negate();
    coupleLower_.negate();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


template<class Type>
void Foam::BlockLduMatrix<Type>::operator=(const BlockLduMatrix<Type>& A)
{
    if (this == &A)
    {
        FatalErrorIn
        (
            "void BlockLduMatrix<Type>::operator="
            "(const BlockLduMatrix<Type>& A)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }

    if (A.lowerPtr_)
    {
        lower() = A.lower();
    }
    else if (lowerPtr_)
    {
        delete lowerPtr_;
        lowerPtr_ = NULL;
    }

    if (A.upperPtr_)
    {
        upper() = A.upper();
    }
    else if (upperPtr_)
    {
        delete upperPtr_;
        upperPtr_ = NULL;
    }

    if (A.diagPtr_)
    {
        diag() = A.diag();
    }

    // Copy interface data
    interfaces_ = A.interfaces_;
    coupleUpper_ = A.coupleUpper_;
    coupleLower_ = A.coupleLower_;

    // Copy constraints
    fixedEqns_ = A.fixedEqns_;
}


template<class Type>
void Foam::BlockLduMatrix<Type>::operator+=(const BlockLduMatrix<Type>& A)
{
    // Do diagonal first
    if (A.thereIsDiag())
    {
        diag() += A.diag();
    }

    if (A.symmetric())
    {
        upper() += A.upper();

        if (this->asymmetric())
        {
            lower() += A.upper().transpose();
        }
    }

    if (A.asymmetric())
    {
        upper() += A.upper();
        lower() += A.lower();
    }

    // Interface data
    coupleUpper_ += A.coupleUpper_;
    coupleLower_ += A.coupleLower_;
}


template<class Type>
void Foam::BlockLduMatrix<Type>::operator-=(const BlockLduMatrix<Type>& A)
{
    // Do diagonal first
    if (A.thereIsDiag())
    {
        diag() -= A.diag();
    }

    if (A.symmetric())
    {
        upper() -= A.upper();

        if (this->asymmetric())
        {
            lower() -= A.upper().transpose();
        }
    }

    if (A.asymmetric())
    {
        upper() -= A.upper();
        lower() -= A.lower();
    }

    // Interface data
    coupleUpper_ -= A.coupleUpper_;
    coupleLower_ -= A.coupleLower_;
}


template<class Type>
void Foam::BlockLduMatrix<Type>::operator*=(const scalarField& sf)
{
    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    //HJ Missing code: add coupling coefficients op*=
    //   HJ, 21/Feb/2008
    // IC - Complicated because we have to scale across the interfaces
    // We may need to include this functionality in lduInterfaceField
    // as additional initInterfaceScale and scaleInterface functions

    if (diagPtr_)
    {
        *diagPtr_ *= sf;
    }

    if (upperPtr_)
    {
        TypeCoeffField& Upper = *upperPtr_;

        const unallocLabelList& l = lduAddr().lowerAddr();

        if (Upper.activeType() == blockCoeffBase::SCALAR)
        {
            scalarTypeField& activeUpper = Upper.asScalar();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeUpper[coeffI] *= sf[l[coeffI]];
            }
        }
        else if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
            linearTypeField& activeUpper = Upper.asLinear();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeUpper[coeffI] *= sf[l[coeffI]];
            }
        }
        else if (Upper.activeType() == blockCoeffBase::SQUARE)
        {
            squareTypeField& activeUpper = Upper.asSquare();

            for (register label coeffI = 0; coeffI < l.size(); coeffI++)
            {
                activeUpper[coeffI] *= sf[l[coeffI]];
            }
        }
    }

    if (lowerPtr_)
    {
        TypeCoeffField& Lower = *lowerPtr_;

        const unallocLabelList& u = lduAddr().upperAddr();

        if (Lower.activeType() == blockCoeffBase::SCALAR)
        {
            scalarTypeField& activeLower = Lower.asScalar();

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                activeLower[coeffI] *= sf[u[coeffI]];
            }
        }
        else if (Lower.activeType() == blockCoeffBase::LINEAR)
        {
            linearTypeField& activeLower = Lower.asLinear();

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                activeLower[coeffI] *= sf[u[coeffI]];
            }
        }
        else if (Lower.activeType() == blockCoeffBase::SQUARE)
        {
            squareTypeField& activeLower = Lower.asSquare();

            for (register label coeffI = 0; coeffI < u.size(); coeffI++)
            {
                activeLower[coeffI] *= sf[u[coeffI]];
            }
        }
    }
}


template<class Type>
void Foam::BlockLduMatrix<Type>::operator*=(const scalar s)
{
    if (diagPtr_)
    {
        *diagPtr_ *= s;
    }

    if (upperPtr_)
    {
        *upperPtr_ *= s;
    }

    if (lowerPtr_)
    {
        *lowerPtr_ *= s;
    }

    // Interface data
    coupleUpper_ *= s;
    coupleLower_ *= s;
}


// ************************************************************************* //
