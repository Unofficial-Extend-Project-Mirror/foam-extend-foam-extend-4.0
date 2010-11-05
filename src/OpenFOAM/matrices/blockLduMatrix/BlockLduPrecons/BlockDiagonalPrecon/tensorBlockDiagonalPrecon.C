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
