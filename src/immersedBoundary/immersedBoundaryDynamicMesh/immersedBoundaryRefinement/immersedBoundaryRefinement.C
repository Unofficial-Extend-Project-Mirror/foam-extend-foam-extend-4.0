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
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryRefinement.H"
#include "fvMesh.H"
#include "immersedBoundaryPolyPatch.H"
#include "triSurfaceDistance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(immersedBoundaryRefinement, 0);
addToRunTimeSelectionTable
(
    refinementSelection,
    immersedBoundaryRefinement,
    dictionary
);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryRefinement::immersedBoundaryRefinement
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    refinementSelection(mesh, dict),
    refinementDistance_
    (
        readScalar(coeffDict().lookup("refinementDistance"))
    ),
    unrefinementDistance_
    (
        readScalar(coeffDict().lookup("unrefinementDistance"))
    ),
    internalRefinementDistance_
    (
        readScalar(coeffDict().lookup("internalRefinementDistance"))
    ),
    internalUnrefinementDistance_
    (
        readScalar(coeffDict().lookup("internalUnrefinementDistance"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::immersedBoundaryRefinement::~immersedBoundaryRefinement()
{}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

Foam::Xfer<Foam::labelList>
Foam::immersedBoundaryRefinement::refinementCellCandidates() const
{
    // Calculate distance to immersed boundary
    // Note: minimum distance is inside (WET) and implies a live cell
    scalarField cellDistance(mesh().nCells(), -GREAT);

    const polyBoundaryMesh& bMesh = mesh().boundaryMesh();

    const vector span(GREAT, GREAT, GREAT);

    forAll (bMesh, patchI)
    {
        if (isA<immersedBoundaryPolyPatch>(bMesh[patchI]))
        {
            const immersedBoundaryPolyPatch& ibPatch =
                refCast<const immersedBoundaryPolyPatch>
                (
                    bMesh[patchI]
                );

            // Create a distance function
            triSurfaceDistance tsd
            (
                ibPatch.triSurfSearch(),
                span,
                ibPatch.internalFlow(),
                false                     // Do not iterate
            );

            // Calculate distance
            cellDistance = Foam::max
            (
                cellDistance,
                tsd.distance(mesh().cellCentres())
            );
        }
    }

    Info<< "Cell distance (min, max): (" << min(cellDistance)
        << ", " << max(cellDistance) << ")" << endl;

    // Create storage for collection of cells. Assume that almost all of the
    // cells will be marked to prevent excessive resizing.
    dynamicLabelList refinementCandidates(mesh().nCells());

    // Create a sorting criterion
    forAll (cellDistance, cellI)
    {
        if
        (
            cellDistance[cellI] > -refinementDistance_
         && cellDistance[cellI] < internalRefinementDistance_
        )
        {
            // Found a refinement cell
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
Foam::immersedBoundaryRefinement::unrefinementPointCandidates() const
{
    // All points are unrefinement candidates (not just split points, these are
    // filtered out in polyhedralRefinement engine)
    dynamicLabelList unrefinementCandidates(mesh().nPoints());

    // Calculate distance to immersed boundary
    // Note: minimum distance is inside (WET) and implies a live cell
    scalarField pointDistance(mesh().nPoints(), -GREAT);

    const polyBoundaryMesh& bMesh = mesh().boundaryMesh();

    const vector span(GREAT, GREAT, GREAT);

    forAll (bMesh, patchI)
    {
        if (isA<immersedBoundaryPolyPatch>(bMesh[patchI]))
        {
            const immersedBoundaryPolyPatch& ibPatch =
                refCast<const immersedBoundaryPolyPatch>(bMesh[patchI]);

            // Create a distance function
            triSurfaceDistance tsd
            (
                ibPatch.triSurfSearch(),
                span,
                ibPatch.internalFlow(),
                false                     // Do not iterate
            );

            // Calculate distance
            pointDistance = Foam::max
            (
                pointDistance,
                tsd.distance(mesh().points())
            );
        }
    }

    Info<< "Split point distance(min, max): (" << min(pointDistance)
        << ", " << max(pointDistance) << ")" << endl;

    forAll (pointDistance, pointI)
    {
        // Unrefine both live cells and dead cells far from the immersed
        // boundary or dead cells
        if
        (
            pointDistance[pointI] < -unrefinementDistance_
         || pointDistance[pointI] > internalUnrefinementDistance_
        )
        {
            unrefinementCandidates.append(pointI);
        }
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
