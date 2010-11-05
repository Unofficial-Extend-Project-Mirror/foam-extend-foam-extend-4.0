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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "BlockCholeskyPrecon.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::BlockCholeskyPrecon<Type>::calcDecoupledPreconDiag()
{
    typedef CoeffField<Type> TypeCoeffField;

    // Note: Assuming lower and upper triangle have the same active type

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
        }
    }
}


template<class Type>
void Foam::BlockCholeskyPrecon<Type>::decoupledPrecondition
(
    Field<Type>& x,
    const Field<Type>& b
) const
{
    typedef CoeffField<Type> TypeCoeffField;

    // Note: Assuming lower and upper triangle have the same active type

    if (this->matrix_.symmetric())
    {
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
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
        }
        else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
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
        }
        else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
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
        }
    }
}


template<class Type>
void Foam::BlockCholeskyPrecon<Type>::decoupledPreconditionT
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
        const TypeCoeffField& LowerCoeff = this->matrix_.lower();
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
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
        }
        else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
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
        }
    }
}


// ************************************************************************* //
