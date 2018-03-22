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

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "minCellSizeRefinement.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(minCellSizeRefinement, 0);
addToRunTimeSelectionTable
(
    refinementSelection,
    minCellSizeRefinement,
    dictionary
);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::minCellSizeRefinement::minCellSizeRefinement
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    refinementSelection(mesh, dict),
    minDelta_(readScalar(coeffDict().lookup("minCellSize")))
{}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::minCellSizeRefinement::~minCellSizeRefinement()
{}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

Foam::Xfer<Foam::labelList>
Foam::minCellSizeRefinement::refinementCellCandidates() const
{
    // Get cell sizes: cube root of cell volume (assuming cube cells)
    const scalarField cellSize = pow(mesh().V().field(), 1.0/3.0);

    // Create storage for collection of cells. Assume that almost all of the
    // cells will be marked to prevent excessive resizing.
    dynamicLabelList refinementCandidates(mesh().nCells());

    // Loop through cells and collect refinement candidates
    forAll (cellSize, cellI)
    {
        if (cellSize[cellI] > minDelta_)
        {
            // Cell is larger than the specified minimum, append cell for
            // potential refinement
            refinementCandidates.append(cellI);
        }
    }

    // Print out some information
    Info<< "Selection algorithm " << type() << " selected "
        << returnReduce(refinementCandidates.size(), sumOp<label>())
        << " cells as refinement candidates."
        << endl;

    // Return the list in the Xfer container to prevent copying
    return refinementCandidates.xfer();
}


Foam::Xfer<Foam::labelList>
Foam::minCellSizeRefinement::unrefinementPointCandidates() const
{
    // Mark all points as unrefinement candidates since only split points may be
    // considered for actual unrefinement and since this refinement criterion
    // will be usually used in combination with others. VV, 15/Mar/2018.

    // All points are unrefinement candidates
    labelList unrefinementCandidates(mesh().nPoints());

    forAll (unrefinementCandidates, pointI)
    {
        unrefinementCandidates[pointI] = pointI;
    }

    // Print out some information
    Info<< "Selection algorithm " << type() << " selected "
        << returnReduce(unrefinementCandidates.size(), sumOp<label>())
        << " points as unrefinement candidates."
        << endl;

    // Return the list in the Xfer container to prevent copying
    return unrefinementCandidates.xfer();
}


// ************************************************************************* //
