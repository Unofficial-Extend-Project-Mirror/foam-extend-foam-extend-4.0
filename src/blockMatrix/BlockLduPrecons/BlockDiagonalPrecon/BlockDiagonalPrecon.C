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
#include "BlockDiagonalPrecon.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::BlockDiagonalPrecon<Type>::precondition
(
    Field<Type>& x,
    const Field<Type>& b
) const
{
    typedef CoeffField<Type> TypeCoeffField;

    typedef typename TypeCoeffField::scalarTypeField scalarTypeField;
    typedef typename TypeCoeffField::linearTypeField linearTypeField;
    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    const TypeCoeffField& diag = this->matrix_.diag();

    if (diag.activeType() == blockCoeffBase::SCALAR)
    {
        const scalarTypeField& activeDiag = diag.asScalar();

        forAll (x, i)
        {
            x[i] = b[i]/activeDiag[i];
        }
    }
    else if (diag.activeType() == blockCoeffBase::LINEAR)
    {
        const linearTypeField& activeDiag = diag.asLinear();

        forAll (x, i)
        {
            x[i] = cmptDivide(b[i], activeDiag[i]);
        }
    }
    else if (diag.activeType() == blockCoeffBase::SQUARE)
    {
        const squareTypeField& activeDiag = diag.asSquare();

        forAll (x, i)
        {
            x[i] = (b[i] & inv(activeDiag[i]));
        }
    }
    else
    {
        FatalErrorIn
        (
            "void BlockDiagonalPrecon<Type>:solve:\n"
            "(\n"
            "    Field<Type>& x,\n"
            "    const Field<Type>& b\n"
            ") const"
        )   << "Problem with coefficient type morphing."
            << abort(FatalError);
    }
}


// ************************************************************************* //
