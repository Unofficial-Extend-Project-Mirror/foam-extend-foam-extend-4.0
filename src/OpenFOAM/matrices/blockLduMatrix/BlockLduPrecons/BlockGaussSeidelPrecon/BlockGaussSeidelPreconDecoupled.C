/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description
    Gauss-Seidel sweep as a preconditioner.  Decoupled version, used in template
    specialisation.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "BlockGaussSeidelPrecon.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

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
        TypeCoeffField dDCoeff = inv(this->matrix_.diag());

        multiply(x, dDCoeff, b);
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
