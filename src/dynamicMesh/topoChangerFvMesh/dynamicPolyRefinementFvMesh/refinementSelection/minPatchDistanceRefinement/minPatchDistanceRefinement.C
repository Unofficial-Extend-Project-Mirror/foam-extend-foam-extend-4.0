/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "minPatchDistanceRefinement.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "patchWave.H"
#include "polyPatchID.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(minPatchDistanceRefinement, 0);
addToRunTimeSelectionTable
(
    refinementSelection,
    minPatchDistanceRefinement,
    dictionary
);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::minPatchDistanceRefinement::minPatchDistanceRefinement
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    refinementSelection(mesh, dict),
    minPatchDistance_(readScalar(coeffDict().lookup("minPatchDistance")))
{}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::minPatchDistanceRefinement::~minPatchDistanceRefinement()
{}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

Foam::Xfer<Foam::labelList>
Foam::minPatchDistanceRefinement::refinementCellCandidates() const
{
    // Note: recalculate distance every time due to probable topo changes

    // Read patch names from dictionary
    const wordList patchNames(coeffDict().lookup("distancePatches"));

    // Collect patchIDs into a hash set
    labelHashSet patchIDs(patchNames.size());
    forAll (patchNames, patchI)
    {
        // Get polyPatchID
        const polyPatchID pID(patchNames[patchI], mesh().boundaryMesh());

        if (pID.active())
        {
            patchIDs.insert(pID.index());
        }
        else
        {
            FatalIOErrorIn
            (
                "Xfer<labelList> minPatchDistanceRefinement::"
                "refinenementPatchCandidates() const",
                coeffDict()
            )   << "Cannot find patch " << patchNames[patchI]
                << " in the mesh." << nl
                << "Available patches are: " << mesh().boundaryMesh().names()
                << abort(FatalIOError);
        }
    }

    // Calculate distance from specified patches and do not correct for accurate
    // near wall distance (false argument)
    patchWave waveDistance(mesh(), patchIDs, false);
    const scalarField& patchDistance = waveDistance.distance();

    // Create storage for collection of cells. Assume that almost all of the
    // cells will be marked to prevent excessive resizing.
    dynamicLabelList refinementCandidates(mesh().nCells());

    // Loop through cells and collect refinement candidates
    forAll (patchDistance, cellI)
    {
        if (patchDistance[cellI] > minPatchDistance_)
        {
            // Cell is far from the patches, append it for potential refinement
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
Foam::minPatchDistanceRefinement::unrefinementPointCandidates() const
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
