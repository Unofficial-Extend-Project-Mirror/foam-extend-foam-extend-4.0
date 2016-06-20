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
