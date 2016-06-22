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

Class
    BlockDiagonalPrecon

Description
    Template specialisation for tensor block diagonal preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef tensorBlockDiagonalPrecon_H
#define tensorBlockDiagonalPrecon_H

#include "BlockDiagonalPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void Foam::BlockDiagonalPrecon<tensor>::precondition
(
    tensorField& x,
    const tensorField& b
) const
{
    typedef CoeffField<tensor> tensorCoeffField;

    typedef tensorCoeffField::scalarTypeField scalarTypeField;
    typedef tensorCoeffField::linearTypeField linearTypeField;

    const tensorCoeffField& diag = matrix_.diag();

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
    else
    {
        FatalErrorIn
        (
            "void BlockDiagonalPrecon<tensor>::solve\n"
            "(\n"
            "    tensorField& x,\n"
            "    const tensorField& b\n"
            ") const"
        )   << "Problem with coefficient type morphing."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
