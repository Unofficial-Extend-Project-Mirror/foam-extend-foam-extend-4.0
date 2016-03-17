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
    Diagonally-corrected Gauss-Seidel sweep as a preconditioner.
    Decoupled version, used in template specialisation.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "BlockDiagGaussSeidelPrecon.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::BlockDiagGaussSeidelPrecon<Type>::calcDecoupledInvDiag()
{
    typedef CoeffField<Type> TypeCoeffField;
    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
    typedef typename TypeCoeffField::scalarType scalarType;

    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::linearType linearType;

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

        if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
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

        if (Upper.activeType() == blockCoeffBase::LINEAR)
        {
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
    else
    {
        FatalErrorIn
        (
            "void BlockDiagGaussSeidelPrecon<Type>::calcDecoupledInvDiag()"
        )   << "Problem with coefficient type morphing."
            << abort(FatalError);
    }
}


template<class Type>
void Foam::BlockDiagGaussSeidelPrecon<Type>::decoupledPrecondition
(
    Field<Type>& x,
    const Field<Type>& b
) const
{
    typedef DecoupledCoeffField<Type> TypeCoeffField;

    if (this->matrix_.diagonal())
    {
        // For a diagonal matrix, the LUDiag_ correction is zero
        // HJ, 26/Jan/2016
        if (invDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            x = invDiag_.asScalar()*b;
        }
        else if (invDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            x = cmptMultiply(invDiag_.asLinear(), b);
        }
    }
    else if (this->matrix_.symmetric() || this->matrix_.asymmetric())
    {
        // Prepare correction
        if (bPlusLU_.empty())
        {
            bPlusLU_.setSize(b.size(), pTraits<Type>::zero);
        }

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

        if (invDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    bPlusLU_ = LUDiag_.asScalar()*x;
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
                    bPlusLU_ = LUDiag_.asScalar()*x;
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
        }
        else if (invDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    bPlusLU_ = cmptMultiply(LUDiag_.asLinear(), x);
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
                    bPlusLU_ = cmptMultiply(LUDiag_.asLinear(), x);
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
        }
        else
        {
            FatalErrorIn
            (
                "void BlockDiagGaussSeidelPrecon<Type>::"
                "decoupledPrecondition\n"
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
            "void BlockDiagGaussSeidelPrecon<Type>::decoupledPrecondition\n"
            "(\n"
            "    Field<Type>& x,\n"
            "    const Field<Type>& b\n"
            ") const"
         )  << "cannot solve incomplete matrix, no diagonal"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::BlockDiagGaussSeidelPrecon<Type>::decoupledPreconditionT
(
    Field<Type>& xT,
    const Field<Type>& bT
) const
{
    typedef DecoupledCoeffField<Type> TypeCoeffField;

    // Prepare correction
    if (bPlusLU_.empty())
    {
        bPlusLU_.setSize(bT.size(), pTraits<Type>::zero);
    }

    if (this->matrix_.diagonal())
    {
        if (invDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            xT = invDiag_.asScalar()*bPlusLU_;
        }
        else if (invDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            xT = cmptMultiply(invDiag_.asLinear(), bPlusLU_);
        }
    }
    else if (this->matrix_.symmetric() || this->matrix_.asymmetric())
    {
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

        if (invDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    bPlusLU_ = LUDiag_.asScalar()*xT;
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
                    bPlusLU_ = LUDiag_.asScalar()*xT;
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
        }
        else if (invDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                for (label sweep = 0; sweep < nSweeps_; sweep++)
                {
                    bPlusLU_ = cmptMultiply(LUDiag_.asLinear(), xT);
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
                    bPlusLU_ = cmptMultiply(LUDiag_.asLinear(), xT);
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
        }
        else
        {
            FatalErrorIn
            (
                "void BlockDiagGaussSeidelPrecon<Type>::"
                "decoupledPreconditionT\n"
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
            "void BlockDiagGaussSeidelPrecon<Type>::decoupledPreconditionT\n"
            "(\n"
            "    Field<Type>& xT,\n"
            "    const Field<Type>& bT\n"
            ") const"
         )  << "cannot solve incomplete matrix, no diagonal"
            << abort(FatalError);
    }
}


// ************************************************************************* //
