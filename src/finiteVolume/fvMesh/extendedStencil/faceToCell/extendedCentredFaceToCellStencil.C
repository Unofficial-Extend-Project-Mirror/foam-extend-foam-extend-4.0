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

#include "mapDistribute.H"
#include "extendedCentredFaceToCellStencil.H"
#include "faceToCellStencil.H"

// Only for access to calcDistributeMap <- needs to be moved out
#include "extendedCellToFaceStencil.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedCentredFaceToCellStencil::extendedCentredFaceToCellStencil
(
    const faceToCellStencil& stencil
)
:
    extendedFaceToCellStencil(stencil.mesh())
{
    stencil_ = stencil;

    // Calculate distribute map (also renumbers elements in stencil)
    mapPtr_ = extendedCellToFaceStencil::calcDistributeMap
    (
        stencil.mesh(),
        stencil.globalNumbering(),
        stencil_
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Per face which elements of the stencil to keep.
void Foam::extendedCentredFaceToCellStencil::compact()
{
    boolList isInStencil(map().constructSize(), false);

    forAll(stencil_, faceI)
    {
        const labelList& stencilCells = stencil_[faceI];

        forAll(stencilCells, i)
        {
            isInStencil[stencilCells[i]] = true;
        }
    }

    mapPtr_().compact(isInStencil);
}


// ************************************************************************* //
