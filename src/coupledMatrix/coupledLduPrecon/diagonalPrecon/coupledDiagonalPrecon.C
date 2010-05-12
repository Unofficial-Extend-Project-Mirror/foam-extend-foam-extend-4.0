/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    coupledDiagonalPrecon

Description
    Diagonal preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*----------------------------------------------------------------------------*/

#include "coupledDiagonalPrecon.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledDiagonalPrecon, 0);

    addToRunTimeSelectionTable
    (
        coupledLduPrecon,
        coupledDiagonalPrecon,
        dictionary
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledDiagonalPrecon::precondition
(
    FieldField<Field, scalar>& x,
    const FieldField<Field, scalar>& b,
    const direction cmpt
) const
{
    // Diagonally precondition all matrices
    forAll (matrix_, rowI)
    {
        scalarField& rowX = x[rowI];
        const scalarField& rowB = b[rowI];
        const scalarField& rowDiag = matrix_[rowI].diag();

        forAll(rowX, i)
        {
            rowX[i] = rowB[i]/rowDiag[i];
        }
    }
}


// ************************************************************************* //
