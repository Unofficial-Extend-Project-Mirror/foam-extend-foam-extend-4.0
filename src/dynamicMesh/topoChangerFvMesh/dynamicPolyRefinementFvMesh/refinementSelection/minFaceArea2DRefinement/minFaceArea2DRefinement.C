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

#include "minFaceArea2DRefinement.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(minFaceArea2DRefinement, 0);
addToRunTimeSelectionTable
(
    refinementSelection,
    minFaceArea2DRefinement,
    dictionary
);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::minFaceArea2DRefinement::minFaceArea2DRefinement
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    refinementSelection(mesh, dict),
    minFaceA_(readScalar(coeffDict().lookup("minFaceArea")))
{
    // It is an error to use this selection algorithm for 3D or 1D cases
    if (mesh.nGeometricD() != 2)
    {
        FatalErrorIn
        (
            "minFaceArea2DRefinement::minFaceArea2DRefinement"
            "\n("
            "\n    const fvMesh& mesh,"
            "\n    const dictionary& dict,"
            "\n)"
        )   << "You are trying to use minFaceArea2DRefinement selection"
            << " algorithm for a case which is not 2D case."
            << nl
            << "This is not allowed."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::minFaceArea2DRefinement::~minFaceArea2DRefinement()
{}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

Foam::Xfer<Foam::labelList>
Foam::minFaceArea2DRefinement::refinementCellCandidates() const
{
    // Get face areas (from polyMesh)
    const scalarField faceArea = mag(mesh().faceAreas());

    // Mark empty or wedge faces
    boolList faceOnEmptyPatch(mesh().nFaces(), false);

    // Loop through all patches
    const polyBoundaryMesh& bMesh = mesh().boundaryMesh();
    forAll (bMesh, patchI)
    {
        const polyPatch& curPatch = bMesh[patchI];

        if (isA<emptyPolyPatch>(curPatch) || isA<wedgePolyPatch>(curPatch))
        {
            // Faces on this patch need to be marked
            const label startFaceI = curPatch.start();
            const label endFaceI = startFaceI + curPatch.size();

            // Loop through all the faces and marke them
            for (label faceI = startFaceI; faceI < endFaceI; ++faceI)
            {
                faceOnEmptyPatch[faceI] = true;
            }
        }
    }


    // Create field that contains maximum face area on empty patch per cell
    scalarField maxEmptyPatchFaceArea(mesh().nCells(), 0.0);

    // Get cells from the mesh
    const cellList& cells = mesh().cells();

    forAll (cells, cellI)
    {
        // Get current cell
        const cell& cellFaces = cells[cellI];

        // Loop through faces of the cell
        forAll (cellFaces, i)
        {
            // Get face index
            const label& faceI = cellFaces[i];

            if (faceOnEmptyPatch[faceI])
            {
                // This face is on empty patch, set maximum face area to the new
                // value if it's larger than the old value
                maxEmptyPatchFaceArea[cellI] =
                    max
                    (
                        maxEmptyPatchFaceArea[cellI],
                        faceArea[faceI]
                    );
            }
        }
    }


    // Create storage for collection of cells. Assume that almost all of the
    // cells will be marked to prevent excessive resizing.
    dynamicLabelList refinementCandidates(mesh().nCells());

    // Loop through cells and collect refinement candidates
    forAll (maxEmptyPatchFaceArea, cellI)
    {
        if (maxEmptyPatchFaceArea[cellI] > minFaceA_)
        {
            // Face area on empty patch is larger than the specified minimum,
            // append cell for potential refinement
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
Foam::minFaceArea2DRefinement::unrefinementPointCandidates() const
{
    // Mark all points as unrefinement candidates since only split points may be
    // considered for actual unrefinement and since this refinement criterion
    // will be usually used in combination with others. VV, 4/Sep/2018.

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
