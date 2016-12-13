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
    Block variant of ILU preconditioning with arbitrary level of fill in (p),
    based on Crout algorithm. Decoupled version.

    Reference: Saad, Y.: Iterative Methods for Sparse Linear Systems
    (2nd Edition), SIAM, 2003.

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "BlockILUCpPrecon.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::BlockILUCpPrecon<Type>::decoupledPrecondition
(
    Field<Type>& x,
    const Field<Type>& b
) const
{
    // Helper type definition of DecoupledCoeffField
    typedef DecoupledCoeffField<Type> DecoupledTypeCoeffField;

    if (!this->matrix_.diagonal())
    {
        // Get upper and lower matrix factors
        const DecoupledTypeCoeffField& Lower = extBlockMatrix_.extendedLower();
        const DecoupledTypeCoeffField& Upper = extBlockMatrix_.extendedUpper();

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
        }
        else
        {
            FatalErrorIn
            (
                "void BlockILUCpPrecon<Type>::decoupledPrecondition\n"
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
            "void BlockILUCpPrecon<Type>::decoupledPrecondition\n"
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
void Foam::BlockILUCpPrecon<Type>::decoupledPreconditionT
(
    Field<Type>& xT,
    const Field<Type>& bT
) const
{
    // Helper type definition of DecoupledCoeffField
    typedef DecoupledCoeffField<Type> DecoupledTypeCoeffField;

    if (!this->matrix_.diagonal())
    {
        // Get upper and lower matrix factors
        const DecoupledTypeCoeffField& Lower = extBlockMatrix_.extendedLower();
        const DecoupledTypeCoeffField& Upper = extBlockMatrix_.extendedUpper();

        // Execute preconditioning by LU substitution.
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
        }
        else
        {
            FatalErrorIn
            (
                "void BlockILUCpPrecon<Type>::decoupledPreconditionT\n"
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
            "void BlockILUCpPrecon<Type>::decoupledPreconditionT\n"
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


// ************************************************************************* //
