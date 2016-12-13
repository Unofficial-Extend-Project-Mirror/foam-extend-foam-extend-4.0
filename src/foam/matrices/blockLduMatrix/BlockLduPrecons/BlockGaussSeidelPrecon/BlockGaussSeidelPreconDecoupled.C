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
    Gauss-Seidel sweep as a preconditioner.  Decoupled version, used in
     template specialisation.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "BlockGaussSeidelPrecon.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::BlockGaussSeidelPrecon<Type>::calcDecoupledInvDiag()
{
    // Get reference to diagonal and obtain inverse by casting
    typedef CoeffField<Type> TypeCoeffField;

    const TypeCoeffField& d = this->matrix_.diag();
    const DecoupledCoeffField<Type>& dd = d;

    invDiag_ = CoeffField<Type>(inv(dd)());
}


template<class Type>
void Foam::BlockGaussSeidelPrecon<Type>::decoupledPrecondition
(
    Field<Type>& x,
    const Field<Type>& b
) const
{
    typedef DecoupledCoeffField<Type> TypeCoeffField;

    if (this->matrix_.diagonal())
    {
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
        }
        else if (invDiag_.activeType() == blockCoeffBase::LINEAR)
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
void Foam::BlockGaussSeidelPrecon<Type>::decoupledPreconditionT
(
    Field<Type>& xT,
    const Field<Type>& bT
) const
{
    typedef DecoupledCoeffField<Type> TypeCoeffField;

    if (this->matrix_.diagonal())
    {
        if (invDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            xT = invDiag_.asScalar()*bT;
        }
        else if (invDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            xT = cmptMultiply(invDiag_.asLinear(), bT);
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
        }
        else if (invDiag_.activeType() == blockCoeffBase::LINEAR)
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
