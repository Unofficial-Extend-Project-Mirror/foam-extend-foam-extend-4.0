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

\*---------------------------------------------------------------------------*/

#include "Stack.H"
#include "triFace.H"
#include "objectMap.H"
#include "changeMap.H"
#include "multiThreader.H"
#include "dynamicTopoFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Method for the bisection of a quad-face in 2D
// - Returns a changeMap with a type specifying:
//     1: Bisection was successful
//    -1: Bisection failed since max number of topo-changes was reached.
//    -2: Bisection failed since resulting quality would be unacceptable.
const changeMap dynamicTopoFvMesh::bisectQuadFace
(
    const label fIndex,
    const changeMap& masterMap,
    bool checkOnly,
    bool forceOp
)
{
    // Quad-face bisection performs the following operations:
    //      [1] Add two points at the middle of the face
    //      [2] Create a new internal face for each bisected cell
    //      [3] Modify existing face and create a new half-face
    //      [4] Modify triangular boundary faces, and create new ones as well
    //      [5] Create edges for new faces
    //      Update faceEdges and edgeFaces information

    // Figure out which thread this is...
    label tIndex = self();

    // Prepare the changeMaps
    changeMap map, slaveMap;

    if
    (
        (statistics_[0] > maxModifications_) &&
        (maxModifications_ > -1)
    )
    {
        // Reached the max allowable topo-changes.
        stack(tIndex).clear();

        return map;
    }

    // Check if edgeRefinements are to be avoided on patch.
    if (lengthEstimator().checkRefinementPatch(whichPatch(fIndex)))
    {
        return map;
    }

    // Sanity check: Is the index legitimate?
    if (fIndex < 0)
    {
        FatalErrorIn
        (
            "\n"
            "const changeMap dynamicTopoFvMesh::bisectQuadFace\n"
            "(\n"
            "    const label fIndex,\n"
            "    const changeMap& masterMap,\n"
            "    bool checkOnly\n"
            ")\n"
        )
            << " Invalid index: " << fIndex << nl
            << " nFaces: " << nFaces_
            << abort(FatalError);
    }

    bool found;
    label replaceFace = -1, retainFace = -1;
    face tmpQuadFace(4), tmpTriFace(3);
    labelList tmpQFEdges(4, -1), tmpTFEdges(3, -1);
    FixedList<label,7> newFaceIndex(-1), newEdgeIndex(-1);
    FixedList<edge,4> commonEdges;
    FixedList<label,4> cornerEdgeIndex(-1), commonEdgeIndex(-1);
    FixedList<label,4> commonFaceIndex(-1);
    FixedList<label,2> newPointIndex(-1), newCellIndex(-1);
    FixedList<label,4> otherEdgeIndex(-1), otherEdgePoint(-1);
    FixedList<label,4> otherPointIndex(-1), nextToOtherPoint(-1);
    FixedList<label,2> c0BdyIndex, c0IntIndex, c1BdyIndex, c1IntIndex;
    FixedList<face,2>  c0BdyFace,  c0IntFace,  c1BdyFace,  c1IntFace;

    // Get the two cells on either side...
    label c0 = owner_[fIndex], c1 = neighbour_[fIndex];

    // Keep track of old / new cells
    FixedList<cell, 2> oldCells(cell(5));
    FixedList<cell, 2> newCells(cell(5));

    // Find the prism faces for cell[0].
    oldCells[0] = cells_[c0];

    meshOps::findPrismFaces
    (
        fIndex,
        c0,
        faces_,
        cells_,
        neighbour_,
        c0BdyFace,
        c0BdyIndex,
        c0IntFace,
        c0IntIndex
    );

    // Check for resulting quality
    if (checkBisection(fIndex, c0BdyIndex[0], forceOp))
    {
        map.type() = -2;
        return map;
    }

    if (c1 != -1)
    {
        // Find the prism faces for cell[1].
        meshOps::findPrismFaces
        (
            fIndex,
            c1,
            faces_,
            cells_,
            neighbour_,
            c1BdyFace,
            c1BdyIndex,
            c1IntFace,
            c1IntIndex
        );

        // Check for resulting quality
        if (checkBisection(fIndex, c1BdyIndex[0], forceOp))
        {
            map.type() = -2;
            return map;
        }
    }

    if (debug > 1)
    {
        Info << nl << nl << "Face: " << fIndex
             << ": " << faces_[fIndex] << " is to be bisected. " << endl;

        label epIndex = whichPatch(fIndex);

        Info << "Patch: ";

        if (epIndex == -1)
        {
            Info << "Internal" << endl;
        }
        else
        {
            Info << boundaryMesh()[epIndex].name() << endl;
        }

        if (debug > 2)
        {
            Info << "Cell[0]: " << c0 << ": " << oldCells[0] << endl;

            forAll(oldCells[0], faceI)
            {
                const labelList& fE = faceEdges_[oldCells[0][faceI]];

                Info << oldCells[0][faceI] << ": "
                     << faces_[oldCells[0][faceI]]
                     << " fE: " << fE
                     << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Info << '\t' << fE[edgeI]
                         << ": " << edges_[fE[edgeI]]
                         << " eF: " << eF
                         << endl;
                }
            }
        }

        // Write out VTK files prior to change
        if (debug > 3)
        {
            labelList cellHull(2, -1);

            cellHull[0] = owner_[fIndex];
            cellHull[1] = neighbour_[fIndex];

            writeVTK
            (
                Foam::name(fIndex)
              + "_Bisect_0",
                cellHull
            );
        }
    }

    // Find the common-edge between the triangular boundary faces
    // and the face under consideration.
    meshOps::findCommonEdge
    (
        c0BdyIndex[0],
        fIndex,
        faceEdges_,
        commonEdgeIndex[0]
    );

    meshOps::findCommonEdge
    (
        c0BdyIndex[1],
        fIndex,
        faceEdges_,
        commonEdgeIndex[1]
    );

    commonEdges[0] = edges_[commonEdgeIndex[0]];
    commonEdges[1] = edges_[commonEdgeIndex[1]];

    // Are we performing only checks?
    if (checkOnly)
    {
        map.type() = 1;
        return map;
    }

    // Find the isolated point on both boundary faces of cell[0]
    meshOps::findIsolatedPoint
    (
        c0BdyFace[0],
        commonEdges[0],
        otherPointIndex[0],
        nextToOtherPoint[0]
    );

    meshOps::findIsolatedPoint
    (
        c0BdyFace[1],
        commonEdges[1],
        otherPointIndex[1],
        nextToOtherPoint[1]
    );

    // For convenience...
    otherEdgePoint[0] = commonEdges[0].otherVertex(nextToOtherPoint[0]);
    otherEdgePoint[1] = commonEdges[1].otherVertex(nextToOtherPoint[1]);

    labelList mP(2, -1);

    // Set mapping for this point
    mP[0] = commonEdges[0][0];
    mP[1] = commonEdges[0][1];

    // Add two new points to the end of the list
    newPointIndex[0] =
    (
        insertPoint
        (
            0.5 * (points_[mP[0]] + points_[mP[1]]),
            0.5 * (oldPoints_[mP[0]] + oldPoints_[mP[1]]),
            mP
        )
    );

    // Set mapping for this point
    mP[0] = commonEdges[1][0];
    mP[1] = commonEdges[1][1];

    newPointIndex[1] =
    (
        insertPoint
        (
            0.5 * (points_[mP[0]] + points_[mP[1]]),
            0.5 * (oldPoints_[mP[0]] + oldPoints_[mP[1]]),
            mP
        )
    );

    {
        map.addPoint(newPointIndex[0]);
        map.addPoint(newPointIndex[1]);
    }

    // Add a new prism cell to the end of the list.
    // Currently invalid, but will be updated later.
    newCellIndex[0] = insertCell(newCells[0], lengthScale_[c0]);

    // Modify the two existing triangle boundary faces

    // Zeroth boundary face - Owner = c[0] & Neighbour [-1] (unchanged)
    meshOps::replaceLabel
    (
        otherEdgePoint[0],
        newPointIndex[0],
        c0BdyFace[0]
    );

    // First boundary face - Owner = newCell[0], Neighbour = -1
    meshOps::replaceLabel
    (
        otherEdgePoint[1],
        newPointIndex[1],
        c0BdyFace[1]
    );

    // Update faces.
    faces_[c0BdyIndex[0]] = c0BdyFace[0];
    faces_[c0BdyIndex[1]] = c0BdyFace[1];

    owner_[c0BdyIndex[1]] = newCellIndex[0];

    meshOps::replaceLabel
    (
        c0BdyIndex[1],
        -1,
        oldCells[0]
    );

    // Detect edges other than commonEdges
    const labelList& fEdges = faceEdges_[fIndex];

    forAll(fEdges, edgeI)
    {
        if
        (
            fEdges[edgeI] != commonEdgeIndex[0] &&
            fEdges[edgeI] != commonEdgeIndex[1]
        )
        {
            // Obtain a reference to this edge
            const edge& eThis = edges_[fEdges[edgeI]];

            if
            (
                eThis[0] == nextToOtherPoint[0]
             || eThis[1] == nextToOtherPoint[0]
            )
            {
                otherEdgeIndex[0] = fEdges[edgeI];
            }
            else
            {
                otherEdgeIndex[1] = fEdges[edgeI];
            }
        }
    }

    // Modify point-labels on the quad face under consideration
    meshOps::replaceLabel
    (
        otherEdgePoint[0],
        newPointIndex[0],
        faces_[fIndex]
    );

    meshOps::replaceLabel
    (
        nextToOtherPoint[1],
        newPointIndex[1],
        faces_[fIndex]
    );

    // Add this face to the map.
    // Although this face isn't technically 'added', it's
    // required for coupled patch mapping.
    map.addFace(fIndex);

    if (debug > 1)
    {
        Info << "Modified face: " << fIndex
             << ": " << faces_[fIndex] << endl;

        if (debug > 2)
        {
            Info << "Common edges: " << nl
                 << commonEdgeIndex[0] << ": " << commonEdges[0] << nl
                 << commonEdgeIndex[1] << ": " << commonEdges[1]
                 << endl;
        }
    }

    // Find the quad face that contains otherEdgeIndex[1]
    found = false;

    const labelList& e1 = faceEdges_[c0IntIndex[0]];

    forAll(e1, edgeI)
    {
        if (e1[edgeI] == otherEdgeIndex[1])
        {
            meshOps::replaceLabel
            (
                c0IntIndex[0],
                -1,
                oldCells[0]
            );

            replaceFace = c0IntIndex[0];
            retainFace = c0IntIndex[1];
            found = true;
            break;
        }
    }

    if (!found)
    {
        // The edge was not found before
        meshOps::replaceLabel
        (
            c0IntIndex[1],
            -1,
            oldCells[0]
        );

        replaceFace = c0IntIndex[1];
        retainFace = c0IntIndex[0];
    }

    // Check if face reversal is necessary for the replacement
    if (owner_[replaceFace] == c0)
    {
        if (neighbour_[replaceFace] == -1)
        {
            // Change the owner
            owner_[replaceFace] = newCellIndex[0];
        }
        else
        {
            // This face has to be reversed
            faces_[replaceFace] = faces_[replaceFace].reverseFace();
            owner_[replaceFace] = neighbour_[replaceFace];
            neighbour_[replaceFace] = newCellIndex[0];

            setFlip(replaceFace);
        }
    }
    else
    {
        // Keep owner, but change neighbour
        neighbour_[replaceFace] = newCellIndex[0];
    }

    // Define the faces for the new cell
    newCells[0][0] = c0BdyIndex[1];
    newCells[0][1] = replaceFace;

    // Define the set of new faces and insert...

    // New interior face; Owner = cell[0] & Neighbour = newCell[0]
    tmpQuadFace[0] = otherPointIndex[0];
    tmpQuadFace[1] = newPointIndex[0];
    tmpQuadFace[2] = newPointIndex[1];
    tmpQuadFace[3] = otherPointIndex[1];

    newFaceIndex[0] =
    (
        insertFace
        (
            -1,
            tmpQuadFace,
            c0,
            newCellIndex[0]
        )
    );

    // Add a faceEdges entry as well
    faceEdges_.append(tmpQFEdges);

    // Find the common edge between quad/quad faces...
    meshOps::findCommonEdge
    (
        c0IntIndex[0],
        c0IntIndex[1],
        faceEdges_,
        otherEdgeIndex[2]
    );

    // ... and size-up edgeFaces for the edge
    meshOps::sizeUpList
    (
        newFaceIndex[0],
        edgeFaces_[otherEdgeIndex[2]]
    );

    meshOps::replaceLabel
    (
        -1,
        newFaceIndex[0],
        newCells[0]
    );

    meshOps::replaceLabel
    (
        -1,
        newFaceIndex[0],
        oldCells[0]
    );

    // remove2DSliver requires this face index for removal
    map.addFace(newFaceIndex[0]);

    // Second boundary face; Owner = newCell[0] & Neighbour = [-1]
    tmpTriFace[0] = otherPointIndex[0];
    tmpTriFace[1] = newPointIndex[0];
    tmpTriFace[2] = otherEdgePoint[0];

    newFaceIndex[1] =
    (
        insertFace
        (
            whichPatch(c0BdyIndex[0]),
            tmpTriFace,
            newCellIndex[0],
            -1
        )
    );

    // Add a faceEdges entry as well
    faceEdges_.append(tmpTFEdges);

    meshOps::replaceLabel
    (
        -1,
        newFaceIndex[1],
        newCells[0]
    );

    // Third boundary face; Owner = c[0] & Neighbour = [-1]
    tmpTriFace[0] = otherPointIndex[1];
    tmpTriFace[1] = newPointIndex[1];
    tmpTriFace[2] = otherEdgePoint[1];

    newFaceIndex[2] =
    (
        insertFace
        (
            whichPatch(c0BdyIndex[1]),
            tmpTriFace,
            c0,
            -1
        )
    );

    // Add a faceEdges entry as well
    faceEdges_.append(tmpTFEdges);

    meshOps::replaceLabel
    (
        -1,
        newFaceIndex[2],
        oldCells[0]
    );

    // Create / modify edges...
    labelList tmpTriEdgeFaces(3, -1);

    // The edge bisecting the zeroth boundary triangular face
    tmpTriEdgeFaces[0] = c0BdyIndex[0];
    tmpTriEdgeFaces[1] = newFaceIndex[0];
    tmpTriEdgeFaces[2] = newFaceIndex[1];

    newEdgeIndex[1] =
    (
        insertEdge
        (
            whichEdgePatch(commonEdgeIndex[0]),
            edge(newPointIndex[0], otherPointIndex[0]),
            tmpTriEdgeFaces
        )
    );

    // Find the common edge between the quad/tri faces...
    meshOps::findCommonEdge
    (
        c0BdyIndex[0],
        replaceFace,
        faceEdges_,
        cornerEdgeIndex[0]
    );

    // ...and correct faceEdges / edgeFaces
    meshOps::replaceLabel
    (
        cornerEdgeIndex[0],
        newEdgeIndex[1],
        faceEdges_[c0BdyIndex[0]]
    );

    meshOps::replaceLabel
    (
        c0BdyIndex[0],
        newFaceIndex[1],
        edgeFaces_[cornerEdgeIndex[0]]
    );

    // The edge bisecting the first boundary triangular face
    tmpTriEdgeFaces[0] = c0BdyIndex[1];
    tmpTriEdgeFaces[1] = newFaceIndex[0];
    tmpTriEdgeFaces[2] = newFaceIndex[2];

    newEdgeIndex[2] =
    (
        insertEdge
        (
            whichEdgePatch(commonEdgeIndex[1]),
            edge(newPointIndex[1], otherPointIndex[1]),
            tmpTriEdgeFaces
        )
    );

    // Find the common edge between the quad/tri faces...
    meshOps::findCommonEdge
    (
        c0BdyIndex[1],
        retainFace,
        faceEdges_,
        cornerEdgeIndex[1]
    );

    // ...and correct faceEdges / edgeFaces
    meshOps::replaceLabel
    (
        cornerEdgeIndex[1],
        newEdgeIndex[2],
        faceEdges_[c0BdyIndex[1]]
    );

    meshOps::replaceLabel
    (
        c0BdyIndex[1],
        newFaceIndex[2],
        edgeFaces_[cornerEdgeIndex[1]]
    );

    if (c1 == -1)
    {
        // The quad boundary face resulting from bisection;
        // Owner = newCell[0] & Neighbour = [-1]
        tmpQuadFace[0] = newPointIndex[1];
        tmpQuadFace[1] = nextToOtherPoint[1];
        tmpQuadFace[2] = otherEdgePoint[0];
        tmpQuadFace[3] = newPointIndex[0];

        newFaceIndex[3] =
        (
            insertFace
            (
                whichPatch(fIndex),
                tmpQuadFace,
                newCellIndex[0],
                -1
            )
        );

        // Add this face to the map.
        map.addFace(newFaceIndex[3]);

        // Add a faceEdges entry as well
        faceEdges_.append(tmpQFEdges);

        // Correct edgeFaces for otherEdgeIndex[1]
        meshOps::replaceLabel
        (
            fIndex,
            newFaceIndex[3],
            edgeFaces_[otherEdgeIndex[1]]
        );

        meshOps::replaceLabel
        (
            -1,
            newFaceIndex[3],
            newCells[0]
        );

        labelList tmpBiEdgeFaces(2, -1);

        // The edge bisecting the face
        tmpTriEdgeFaces[0] = newFaceIndex[3];
        tmpTriEdgeFaces[1] = newFaceIndex[0];
        tmpTriEdgeFaces[2] = fIndex;

        newEdgeIndex[0] =
        (
            insertEdge
            (
                whichEdgePatch(otherEdgeIndex[0]),
                edge(newPointIndex[0], newPointIndex[1]),
                tmpTriEdgeFaces
            )
        );

        // Replace an edge on the bisected face
        meshOps::replaceLabel
        (
            otherEdgeIndex[1],
            newEdgeIndex[0],
            faceEdges_[fIndex]
        );

        // Create / replace side edges created from face bisection
        tmpBiEdgeFaces[0] = newFaceIndex[1];
        tmpBiEdgeFaces[1] = newFaceIndex[3];

        newEdgeIndex[3] =
        (
            insertEdge
            (
                whichEdgePatch(commonEdgeIndex[0]),
                edge(newPointIndex[0], otherEdgePoint[0]),
                tmpBiEdgeFaces
            )
        );

        tmpBiEdgeFaces[0] = c0BdyIndex[1];
        tmpBiEdgeFaces[1] = newFaceIndex[3];

        newEdgeIndex[4] =
        (
            insertEdge
            (
                whichEdgePatch(commonEdgeIndex[1]),
                edge(newPointIndex[1], nextToOtherPoint[1]),
                tmpBiEdgeFaces
            )
        );

        // Now that edges are defined, configure faceEdges
        // for all new faces

        // The quad interior face; Owner = cell[0] & Neighbour = newCell[0]
        tmpQFEdges[0] = newEdgeIndex[0];
        tmpQFEdges[1] = newEdgeIndex[1];
        tmpQFEdges[2] = otherEdgeIndex[2];
        tmpQFEdges[3] = newEdgeIndex[2];
        faceEdges_[newFaceIndex[0]] = tmpQFEdges;

        // Second boundary face; Owner = newCell[0] & Neighbour = [-1]
        tmpTFEdges[0] = newEdgeIndex[3];
        tmpTFEdges[1] = cornerEdgeIndex[0];
        tmpTFEdges[2] = newEdgeIndex[1];
        faceEdges_[newFaceIndex[1]] = tmpTFEdges;

        // Third boundary face; Owner = c[0] & Neighbour = [-1]
        tmpTFEdges[0] = newEdgeIndex[2];
        tmpTFEdges[1] = cornerEdgeIndex[1];
        tmpTFEdges[2] = commonEdgeIndex[1];
        faceEdges_[newFaceIndex[2]] = tmpTFEdges;

        // The quad face from bisection:
        tmpQFEdges[0] = otherEdgeIndex[1];
        tmpQFEdges[1] = newEdgeIndex[3];
        tmpQFEdges[2] = newEdgeIndex[0];
        tmpQFEdges[3] = newEdgeIndex[4];
        faceEdges_[newFaceIndex[3]] = tmpQFEdges;

        meshOps::replaceLabel
        (
            commonEdgeIndex[1],
            newEdgeIndex[4],
            faceEdges_[c0BdyIndex[1]]
        );

        meshOps::replaceLabel
        (
            c0BdyIndex[1],
            newFaceIndex[2],
            edgeFaces_[commonEdgeIndex[1]]
        );

        if (debug > 2)
        {
            Info << "Modified Cell[0]: "
                 << c0 << ": " << oldCells[0] << endl;

            forAll(oldCells[0], faceI)
            {
                const labelList& fE = faceEdges_[oldCells[0][faceI]];

                Info << oldCells[0][faceI]
                     << ": " << faces_[oldCells[0][faceI]]
                     << " fE: " << fE
                     << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Info << '\t' << fE[edgeI]
                         << ": " << edges_[fE[edgeI]]
                         << " eF: " << eF
                         << endl;
                }
            }

            Info << "New Cell[0]: " << newCellIndex[0]
                 << ": " << newCells[0] << endl;

            forAll(newCells[0], faceI)
            {
                const labelList& fE = faceEdges_[newCells[0][faceI]];

                Info << newCells[0][faceI]
                     << ": " << faces_[newCells[0][faceI]]
                     << " fE: " << fE
                     << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Info << '\t' << fE[edgeI]
                         << ": " << edges_[fE[edgeI]]
                         << " eF: " << eF
                         << endl;
                }
            }
        }
    }
    else
    {
        oldCells[1] = cells_[c1];

        newCellIndex[1] = insertCell(newCells[1], lengthScale_[c1]);

        if (debug > 2)
        {
            Info << "Cell[1]: " << c1 << ": " << oldCells[1] << endl;

            forAll(oldCells[1], faceI)
            {
                const labelList& fE = faceEdges_[oldCells[1][faceI]];

                Info << oldCells[1][faceI] << ": "
                     << faces_[oldCells[1][faceI]]
                     << " fE: " << fE
                     << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Info << '\t' << fE[edgeI]
                         << ": " << edges_[fE[edgeI]]
                         << " eF: " << eF
                         << endl;
                }
            }
        }

        // Find the interior face that contains otherEdgeIndex[1]
        found = false;

        const labelList& e2 = faceEdges_[c1IntIndex[0]];

        forAll(e2, edgeI)
        {
            if (e2[edgeI] == otherEdgeIndex[1])
            {
                meshOps::replaceLabel
                (
                    c1IntIndex[0],
                    -1,
                    oldCells[1]
                );

                replaceFace = c1IntIndex[0];
                retainFace = c1IntIndex[1];
                found = true;
                break;
            }
        }

        if (!found)
        {
            // The edge was not found before
            meshOps::replaceLabel
            (
                c1IntIndex[1],
                -1,
                oldCells[1]
            );

            replaceFace = c1IntIndex[1];
            retainFace = c1IntIndex[0];
        }

        // Check if face reversal is necessary for the replacement
        if (owner_[replaceFace] == c1)
        {
            if (neighbour_[replaceFace] == -1)
            {
                // Change the owner
                owner_[replaceFace] = newCellIndex[1];
            }
            else
            {
                // This face has to be reversed
                faces_[replaceFace] = faces_[replaceFace].reverseFace();
                owner_[replaceFace] = neighbour_[replaceFace];
                neighbour_[replaceFace] = newCellIndex[1];

                setFlip(replaceFace);
            }
        }
        else
        {
            // Keep owner, but change neighbour
            neighbour_[replaceFace] = newCellIndex[1];
        }

        // Define attributes for the new prism cell
        newCells[1][0] = replaceFace;

        // The quad interior face resulting from bisection;
        // Owner = newCell[0] & Neighbour = newCell[1]
        tmpQuadFace[0] = newPointIndex[1];
        tmpQuadFace[1] = nextToOtherPoint[1];
        tmpQuadFace[2] = otherEdgePoint[0];
        tmpQuadFace[3] = newPointIndex[0];

        newFaceIndex[3] =
        (
            insertFace
            (
                -1,
                tmpQuadFace,
                newCellIndex[0],
                newCellIndex[1]
            )
        );

        // Add a faceEdges entry as well
        faceEdges_.append(tmpQFEdges);

        // Correct edgeFaces for otherEdgeIndex[1]
        meshOps::replaceLabel
        (
            fIndex,
            newFaceIndex[3],
            edgeFaces_[otherEdgeIndex[1]]
        );

        meshOps::replaceLabel
        (
            -1,
            newFaceIndex[3],
            newCells[0]
        );

        meshOps::replaceLabel
        (
            -1,
            newFaceIndex[3],
            newCells[1]
        );

        newCells[1][1] = newFaceIndex[3];

        // Check for common edges among the two boundary faces
        // Find the isolated point on both boundary faces of cell[1]
        if
        (
            meshOps::findCommonEdge
            (
                c1BdyIndex[0],
                c0BdyIndex[0],
                faceEdges_,
                commonEdgeIndex[2]
            )
        )
        {
            meshOps::findCommonEdge
            (
                c1BdyIndex[1],
                c0BdyIndex[1],
                faceEdges_,
                commonEdgeIndex[3]
            );

            commonFaceIndex[2] = c1BdyIndex[0];
            commonFaceIndex[3] = c1BdyIndex[1];
        }
        else
        {
            meshOps::findCommonEdge
            (
                c1BdyIndex[0],
                c0BdyIndex[1],
                faceEdges_,
                commonEdgeIndex[3]
            );

            meshOps::findCommonEdge
            (
                c1BdyIndex[1],
                c0BdyIndex[0],
                faceEdges_,
                commonEdgeIndex[2]
            );

            commonFaceIndex[2] = c1BdyIndex[1];
            commonFaceIndex[3] = c1BdyIndex[0];
        }

        commonEdges[2] = edges_[commonEdgeIndex[2]];
        commonEdges[3] = edges_[commonEdgeIndex[3]];

        if (debug > 2)
        {
            Info << "Common edges: " << nl
                 << commonEdgeIndex[2] << ": " << commonEdges[2] << nl
                 << commonEdgeIndex[3] << ": " << commonEdges[3]
                 << endl;
        }

        meshOps::findIsolatedPoint
        (
            faces_[commonFaceIndex[2]],
            commonEdges[2],
            otherPointIndex[2],
            nextToOtherPoint[2]
        );

        meshOps::findIsolatedPoint
        (
            faces_[commonFaceIndex[3]],
            commonEdges[3],
            otherPointIndex[3],
            nextToOtherPoint[3]
        );

        // For convenience...
        otherEdgePoint[2] = commonEdges[2].otherVertex(nextToOtherPoint[2]);
        otherEdgePoint[3] = commonEdges[3].otherVertex(nextToOtherPoint[3]);

        // Modify the two existing triangle boundary faces

        // Zeroth boundary face - Owner = newCell[1], Neighbour = -1
        meshOps::replaceLabel
        (
            otherEdgePoint[2],
            newPointIndex[0],
            faces_[commonFaceIndex[2]]
        );

        owner_[commonFaceIndex[2]] = newCellIndex[1];

        meshOps::replaceLabel
        (
            commonFaceIndex[2],
            -1,
            oldCells[1]
        );

        newCells[1][2] = commonFaceIndex[2];

        // First boundary face - Owner = c[1] & Neighbour [-1] (unchanged)
        meshOps::replaceLabel
        (
            otherEdgePoint[3],
            newPointIndex[1],
            faces_[commonFaceIndex[3]]
        );

        // New interior face; Owner = cell[1] & Neighbour = newCell[1]
        tmpQuadFace[0] = newPointIndex[0];
        tmpQuadFace[1] = otherPointIndex[2];
        tmpQuadFace[2] = otherPointIndex[3];
        tmpQuadFace[3] = newPointIndex[1];

        newFaceIndex[4] =
        (
            insertFace
            (
                -1,
                tmpQuadFace,
                c1,
                newCellIndex[1]
            )
        );

        // Add a faceEdges entry as well
        faceEdges_.append(tmpQFEdges);

        // remove2DSliver requires this face index for removal
        map.addFace(newFaceIndex[4]);

        // Find the common edge between quad/quad faces...
        meshOps::findCommonEdge
        (
            c1IntIndex[0],
            c1IntIndex[1],
            faceEdges_,
            otherEdgeIndex[3]
        );

        // ... and size-up edgeFaces for the edge
        meshOps::sizeUpList
        (
            newFaceIndex[4],
            edgeFaces_[otherEdgeIndex[3]]
        );

        meshOps::replaceLabel
        (
            -1,
            newFaceIndex[4],
            newCells[1]
        );

        meshOps::replaceLabel
        (
            -1,
            newFaceIndex[4],
            oldCells[1]
        );

        // Second boundary face; Owner = cell[1] & Neighbour [-1]
        tmpTriFace[0] = otherPointIndex[2];
        tmpTriFace[1] = newPointIndex[0];
        tmpTriFace[2] = otherEdgePoint[2];

        newFaceIndex[5] =
        (
            insertFace
            (
                whichPatch(commonFaceIndex[2]),
                tmpTriFace,
                c1,
                -1
            )
        );

        // Add a faceEdges entry as well
        faceEdges_.append(tmpTFEdges);

        meshOps::replaceLabel
        (
            -1,
            newFaceIndex[5],
            oldCells[1]
        );

        // Third boundary face; Owner = newCell[1] & Neighbour [-1]
        tmpTriFace[0] = otherPointIndex[3];
        tmpTriFace[1] = newPointIndex[1];
        tmpTriFace[2] = otherEdgePoint[3];

        newFaceIndex[6] =
        (
            insertFace
            (
                whichPatch(commonFaceIndex[3]),
                tmpTriFace,
                newCellIndex[1],
                -1
            )
        );

        // Add a faceEdges entry as well
        faceEdges_.append(tmpTFEdges);

        meshOps::replaceLabel
        (
            -1,
            newFaceIndex[6],
            newCells[1]
        );

        // Create / modify edges...
        labelList tmpQuadEdgeFaces(4, -1);

        // The internal edge bisecting the face
        tmpQuadEdgeFaces[0] = fIndex;
        tmpQuadEdgeFaces[1] = newFaceIndex[0];
        tmpQuadEdgeFaces[2] = newFaceIndex[3];
        tmpQuadEdgeFaces[3] = newFaceIndex[4];

        newEdgeIndex[0] =
        (
            insertEdge
            (
                -1,
                edge(newPointIndex[0], newPointIndex[1]),
                tmpQuadEdgeFaces
            )
        );

        // Replace an edge on the bisected face
        meshOps::replaceLabel
        (
            otherEdgeIndex[1],
            newEdgeIndex[0],
            faceEdges_[fIndex]
        );

        // Create / replace side edges created from face bisection
        tmpTriEdgeFaces[0] = commonFaceIndex[2];
        tmpTriEdgeFaces[1] = newFaceIndex[3];
        tmpTriEdgeFaces[2] = newFaceIndex[1];

        newEdgeIndex[3] =
        (
            insertEdge
            (
                whichEdgePatch(commonEdgeIndex[2]),
                edge(newPointIndex[0], nextToOtherPoint[2]),
                tmpTriEdgeFaces
            )
        );

        tmpTriEdgeFaces[0] = c0BdyIndex[1];
        tmpTriEdgeFaces[1] = newFaceIndex[3];
        tmpTriEdgeFaces[2] = newFaceIndex[6];

        newEdgeIndex[4] =
        (
            insertEdge
            (
                whichEdgePatch(commonEdgeIndex[3]),
                edge(newPointIndex[1], otherEdgePoint[3]),
                tmpTriEdgeFaces
            )
        );

        // The edge bisecting the second boundary triangular face
        tmpTriEdgeFaces[0] = commonFaceIndex[2];
        tmpTriEdgeFaces[1] = newFaceIndex[4];
        tmpTriEdgeFaces[2] = newFaceIndex[5];

        newEdgeIndex[5] =
        (
            insertEdge
            (
                whichEdgePatch(commonEdgeIndex[2]),
                edge(newPointIndex[0], otherPointIndex[2]),
                tmpTriEdgeFaces
            )
        );

        // Find the common edge between the quad/tri faces...
        meshOps::findCommonEdge
        (
            commonFaceIndex[2],
            retainFace,
            faceEdges_,
            cornerEdgeIndex[2]
        );

        // ...and correct faceEdges / edgeFaces
        meshOps::replaceLabel
        (
            cornerEdgeIndex[2],
            newEdgeIndex[5],
            faceEdges_[commonFaceIndex[2]]
        );

        meshOps::replaceLabel
        (
            commonFaceIndex[2],
            newFaceIndex[5],
            edgeFaces_[cornerEdgeIndex[2]]
        );

        // The edge bisecting the third boundary triangular face
        tmpTriEdgeFaces[0] = commonFaceIndex[3];
        tmpTriEdgeFaces[1] = newFaceIndex[4];
        tmpTriEdgeFaces[2] = newFaceIndex[6];

        newEdgeIndex[6] =
        (
            insertEdge
            (
                whichEdgePatch(commonEdgeIndex[3]),
                edge(newPointIndex[1], otherPointIndex[3]),
                tmpTriEdgeFaces
            )
        );

        // Find the common edge between the quad/tri faces...
        meshOps::findCommonEdge
        (
            commonFaceIndex[3],
            replaceFace,
            faceEdges_,
            cornerEdgeIndex[3]
        );

        // ...and correct faceEdges / edgeFaces
        meshOps::replaceLabel
        (
            cornerEdgeIndex[3],
            newEdgeIndex[6],
            faceEdges_[commonFaceIndex[3]]
        );

        meshOps::replaceLabel
        (
            commonFaceIndex[3],
            newFaceIndex[6],
            edgeFaces_[cornerEdgeIndex[3]]
        );

        // Now that edges are defined, configure faceEdges
        // for all new faces

        // The quad interior face; Owner = c[0] & Neighbour = newCell[0]
        tmpQFEdges[0] = newEdgeIndex[0];
        tmpQFEdges[1] = newEdgeIndex[1];
        tmpQFEdges[2] = otherEdgeIndex[2];
        tmpQFEdges[3] = newEdgeIndex[2];
        faceEdges_[newFaceIndex[0]] = tmpQFEdges;

        // Second boundary face; Owner = newCell[0] & Neighbour = [-1]
        tmpTFEdges[0] = newEdgeIndex[3];
        tmpTFEdges[1] = cornerEdgeIndex[0];
        tmpTFEdges[2] = newEdgeIndex[1];
        faceEdges_[newFaceIndex[1]] = tmpTFEdges;

        // Third boundary face; Owner = c[0] & Neighbour = [-1]
        tmpTFEdges[0] = newEdgeIndex[2];
        tmpTFEdges[1] = cornerEdgeIndex[1];
        tmpTFEdges[2] = commonEdgeIndex[3];
        faceEdges_[newFaceIndex[2]] = tmpTFEdges;

        // The quad face from bisection:
        tmpQFEdges[0] = otherEdgeIndex[1];
        tmpQFEdges[1] = newEdgeIndex[3];
        tmpQFEdges[2] = newEdgeIndex[0];
        tmpQFEdges[3] = newEdgeIndex[4];
        faceEdges_[newFaceIndex[3]] = tmpQFEdges;

        // The quad interior face; Owner = c[1] & Neighbour = newCell[1]
        tmpQFEdges[0] = newEdgeIndex[0];
        tmpQFEdges[1] = newEdgeIndex[5];
        tmpQFEdges[2] = otherEdgeIndex[3];
        tmpQFEdges[3] = newEdgeIndex[6];
        faceEdges_[newFaceIndex[4]] = tmpQFEdges;

        // Second boundary face; Owner = c[1] & Neighbour = [-1]
        tmpTFEdges[0] = commonEdgeIndex[2];
        tmpTFEdges[1] = cornerEdgeIndex[2];
        tmpTFEdges[2] = newEdgeIndex[5];
        faceEdges_[newFaceIndex[5]] = tmpTFEdges;

        // Third boundary face; Owner = newCell[1] & Neighbour = [-1]
        tmpTFEdges[0] = newEdgeIndex[4];
        tmpTFEdges[1] = cornerEdgeIndex[3];
        tmpTFEdges[2] = newEdgeIndex[6];
        faceEdges_[newFaceIndex[6]] = tmpTFEdges;

        meshOps::replaceLabel
        (
            commonEdgeIndex[1],
            newEdgeIndex[4],
            faceEdges_[c0BdyIndex[1]]
        );

        meshOps::replaceLabel
        (
            c0BdyIndex[1],
            newFaceIndex[2],
            edgeFaces_[commonEdgeIndex[1]]
        );

        meshOps::replaceLabel
        (
            commonEdgeIndex[2],
            newEdgeIndex[3],
            faceEdges_[commonFaceIndex[2]]
        );

        meshOps::replaceLabel
        (
            commonFaceIndex[2],
            newFaceIndex[5],
            edgeFaces_[commonEdgeIndex[2]]
        );

        if (debug > 2)
        {
            Info << nl << "Modified Cell[0]: "
                 << c0 << ": " << oldCells[0] << endl;

            forAll(oldCells[0], faceI)
            {
                const labelList& fE = faceEdges_[oldCells[0][faceI]];

                Info << oldCells[0][faceI]
                     << ": " << faces_[oldCells[0][faceI]]
                     << " fE: " << fE
                     << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Info << '\t' << fE[edgeI]
                         << ": " << edges_[fE[edgeI]]
                         << " eF: " << eF
                         << endl;
                }
            }

            Info << "New Cell[0]: "
                 << newCellIndex[0] << ": " << newCells[0] << endl;

            forAll(newCells[0], faceI)
            {
                const labelList& fE = faceEdges_[newCells[0][faceI]];

                Info << newCells[0][faceI] << ": "
                     << faces_[newCells[0][faceI]]
                     << " fE: " << fE
                     << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Info << '\t' << fE[edgeI]
                         << ": " << edges_[fE[edgeI]]
                         << " eF: " << eF
                         << endl;
                }
            }

            Info << nl << "Modified Cell[1]: "
                 << c1 << ": " << oldCells[1] << endl;

            forAll(oldCells[1], faceI)
            {
                const labelList& fE = faceEdges_[oldCells[1][faceI]];

                Info << oldCells[1][faceI] << ": "
                     << faces_[oldCells[1][faceI]]
                     << " fE: " << fE
                     << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Info << '\t' << fE[edgeI]
                         << ": " << edges_[fE[edgeI]]
                         << " eF: " << eF
                         << endl;
                }
            }

            Info << "New Cell[1]: "
                 << newCellIndex[1] << ": " << newCells[1] << endl;

            forAll(newCells[1], faceI)
            {
                const labelList& fE = faceEdges_[newCells[1][faceI]];

                Info << newCells[1][faceI] << ": "
                     << faces_[newCells[1][faceI]]
                     << " fE: " << fE
                     << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Info << '\t' << fE[edgeI]
                         << ": " << edges_[fE[edgeI]]
                         << " eF: " << eF
                         << endl;
                }
            }
        }

        // Update the cell list.
        cells_[c1] = oldCells[1];
        cells_[newCellIndex[1]] = newCells[1];
    }

    // Update the cell list.
    cells_[c0] = oldCells[0];
    cells_[newCellIndex[0]] = newCells[0];

    // Modify point labels for common edges
    if (edges_[commonEdgeIndex[0]].start() == otherEdgePoint[0])
    {
        edges_[commonEdgeIndex[0]].start() = newPointIndex[0];
    }
    else
    {
        edges_[commonEdgeIndex[0]].end() = newPointIndex[0];
    }

    if (edges_[commonEdgeIndex[1]].start() == nextToOtherPoint[1])
    {
        edges_[commonEdgeIndex[1]].start() = newPointIndex[1];
    }
    else
    {
        edges_[commonEdgeIndex[1]].end() = newPointIndex[1];
    }

    // Write out VTK files after change
    if (debug > 3)
    {
        labelList cellHull(4, -1);

        cellHull[0] = owner_[fIndex];
        cellHull[1] = neighbour_[fIndex];
        cellHull[2] = owner_[newFaceIndex[3]];
        cellHull[3] = neighbour_[newFaceIndex[3]];

        writeVTK
        (
            Foam::name(fIndex)
          + "_Bisect_1",
            cellHull
        );
    }

    // Fill-in mapping information
    FixedList<label, 4> mapCells(-1);

    mapCells[0] = c0;
    mapCells[1] = newCellIndex[0];

    if (c1 != -1)
    {
        mapCells[2] = c1;
        mapCells[3] = newCellIndex[1];
    }

    labelList mC(1, c0);

    forAll(mapCells, cellI)
    {
        if (mapCells[cellI] == -1)
        {
            continue;
        }

        // Set the mapping for this cell
        setCellMapping(mapCells[cellI], mC);
    }

    // Set fill-in mapping information for the modified face.
    if (c1 == -1)
    {
        // Set the mapping for this face
        setFaceMapping(fIndex, labelList(1, fIndex));
    }
    else
    {
        // Internal face. Default mapping.
        setFaceMapping(fIndex);
    }

    forAll(newFaceIndex, faceI)
    {
        if (newFaceIndex[faceI] == -1)
        {
            continue;
        }

        // Check for boundary faces
        if (neighbour_[newFaceIndex[faceI]] == -1)
        {
            // Boundary face. Compute mapping.
            labelList mC;

            if (faces_[newFaceIndex[faceI]].size() == 4)
            {
                // Quad-face on boundary
                mC.setSize(1, fIndex);
            }
            else
            if (faces_[newFaceIndex[faceI]].size() == 3)
            {
                label triFacePatch = whichPatch(newFaceIndex[faceI]);

                // Fetch face-normals
                vector tfNorm, f0Norm, f1Norm;

                tfNorm = faces_[newFaceIndex[faceI]].normal(oldPoints_);
                f0Norm = faces_[c0BdyIndex[0]].normal(oldPoints_);
                f1Norm = faces_[c0BdyIndex[1]].normal(oldPoints_);

                // Tri-face on boundary. Perform normal checks
                // also, because of empty patches.
                if
                (
                    (whichPatch(c0BdyIndex[0]) == triFacePatch) &&
                    ((tfNorm & f0Norm) > 0.0)
                )
                {
                    mC.setSize(1, c0BdyIndex[0]);
                }
                else
                if
                (
                    (whichPatch(c0BdyIndex[1]) == triFacePatch) &&
                    ((tfNorm & f1Norm) > 0.0)
                )
                {
                    mC.setSize(1, c0BdyIndex[1]);
                }
                else
                {
                    FatalErrorIn
                    (
                        "\n"
                        "const changeMap dynamicTopoFvMesh::bisectQuadFace\n"
                        "(\n"
                        "    const label fIndex,\n"
                        "    const changeMap& masterMap,\n"
                        "    bool checkOnly\n"
                        ")\n"
                    )
                        << " Unable to find patch for face: "
                        << newFaceIndex[faceI] << ":: "
                        << faces_[newFaceIndex[faceI]] << nl
                        << " Patch: " << triFacePatch << nl
                        << abort(FatalError);
                }
            }

            // Set the mapping for this face
            setFaceMapping(newFaceIndex[faceI], mC);
        }
        else
        {
            // Internal quad-faces get default mapping.
            setFaceMapping(newFaceIndex[faceI]);
        }
    }

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    statistics_[3]++;

    // Increment surface-counter
    if (c1 == -1)
    {
        statistics_[5]++;
    }

    // Increment the number of modifications
    statistics_[0]++;

    // Specify that the operation was successful
    map.type() = 1;

    // Return the changeMap
    return map;
}

// Method for the bisection of an edge in 3D
// - Returns a changeMap with a type specifying:
//     1: Bisection was successful
//    -1: Bisection failed since max number of topo-changes was reached.
//    -2: Bisection failed since resulting quality would be unacceptable.
// - AddedPoints contain the index of the newly added point.
const changeMap dynamicTopoFvMesh::bisectEdge
(
    const label eIndex,
    bool checkOnly,
    bool forceOp
)
{
    // Edge bisection performs the following operations:
    //      [1] Add a point at middle of the edge
    //      [2] Bisect all faces surrounding this edge
    //      [3] Bisect all cells surrounding this edge
    //      [4] Create internal/external edges for each bisected face
    //      [5] Create internal faces for each bisected cell
    //      Update faceEdges, edgeFaces and edgePoints information

    // For 2D meshes, perform face-bisection
    if (twoDMesh_)
    {
        return bisectQuadFace(eIndex, changeMap(), checkOnly);
    }

    // Figure out which thread this is...
    label tIndex = self();

    // Prepare the changeMaps
    changeMap map, slaveMap;

    if
    (
        (statistics_[0] > maxModifications_) &&
        (maxModifications_ > -1)
    )
    {
        // Reached the max allowable topo-changes.
        stack(tIndex).clear();

        return map;
    }

    // Check if edgeRefinements are to be avoided on patch.
    if (lengthEstimator().checkRefinementPatch(whichEdgePatch(eIndex)))
    {
        return map;
    }

    // Sanity check: Is the index legitimate?
    if (eIndex < 0)
    {
        FatalErrorIn
        (
            "\n"
            "const changeMap dynamicTopoFvMesh::bisectEdge\n"
            "(\n"
            "    const label eIndex,\n"
            "    bool checkOnly,\n"
            "    bool forceOp\n"
            ")\n"
        )
            << " Invalid index: " << eIndex
            << " nEdges: " << nEdges_
            << abort(FatalError);
    }

    // Before we bisect this edge, check whether the operation will
    // yield an acceptable cell-quality.
    scalar minQ = 0.0;

    if ((minQ = computeBisectionQuality(eIndex)) < sliverThreshold_)
    {
        // Check if the quality is actually valid before forcing it.
        if (forceOp && (minQ < 0.0))
        {
            FatalErrorIn
            (
                "\n"
                "const changeMap dynamicTopoFvMesh::bisectEdge\n"
                "(\n"
                "    const label eIndex,\n"
                "    bool checkOnly,\n"
                "    bool forceOp\n"
                ")\n"
            )
                << " Forcing bisection on edge: " << eIndex
                << " will yield an invalid cell."
                << abort(FatalError);
        }
        else
        if (!forceOp)
        {
            map.type() = -2;
            return map;
        }
    }

    // Are we performing only checks?
    if (checkOnly)
    {
        map.type() = 1;
        return map;
    }

    // Update number of surface bisections, if necessary.
    if (whichEdgePatch(eIndex) > -1)
    {
        statistics_[5]++;
    }

    // Hull variables
    face tmpTriFace(3);
    labelList tmpEdgeFaces(3,-1);
    labelList tmpIntEdgeFaces(4,-1);
    labelList tmpEdgePoints(3,-1);
    labelList tmpIntEdgePoints(4,-1);
    labelList tmpFaceEdges(3,-1);

    // Make a copy of existing entities
    const labelList vertexHull = edgePoints_[eIndex];
    label m = vertexHull.size();

    // Size up the hull lists
    labelList cellHull(m, -1);
    labelList faceHull(m, -1);
    labelList edgeHull(m, -1);
    labelListList ringEntities(4, labelList(m, -1));

    // Construct a hull around this edge
    meshOps::constructHull
    (
        eIndex,
        faces_,
        edges_,
        cells_,
        owner_,
        neighbour_,
        faceEdges_,
        edgeFaces_,
        edgePoints_,
        edgeHull,
        faceHull,
        cellHull,
        ringEntities
    );

    if (debug > 1)
    {
        Info << nl << nl
             << "Edge: " << eIndex
             << ": " << edges_[eIndex]
             << " is to be bisected. " << endl;

        label epIndex = whichEdgePatch(eIndex);

        Info << "Patch: ";

        if (epIndex == -1)
        {
            Info << "Internal" << endl;
        }
        else
        {
            Info << boundaryMesh()[epIndex].name() << endl;
        }

        // Write out VTK files prior to change
        if (debug > 3)
        {
            writeVTK
            (
                Foam::name(eIndex)
              + "_Bisect_0",
                cellHull
            );
        }
    }

    labelList mP(2, -1);

    // Set mapping for this point
    mP[0] = edges_[eIndex][0];
    mP[1] = edges_[eIndex][1];

    // Add a new point to the end of the list
    label newPointIndex =
    (
        insertPoint
        (
            0.5 * (points_[mP[0]] + points_[mP[1]]),
            0.5 * (oldPoints_[mP[0]] + oldPoints_[mP[1]]),
            mP
        )
    );

    // Add this point to the map.
    map.addPoint(newPointIndex);

    // New edges can lie on a bounding curve between
    // coupled and non-coupled faces. Preferentially
    // add edges to coupled-patches, if they exist.
    // This makes it convenient for coupled patch matching.
    label nePatch = -1;

    {
        nePatch = whichEdgePatch(eIndex);
    }

    // Add a new edge to the end of the list
    label newEdgeIndex =
    (
        insertEdge
        (
            nePatch,
            edge(newPointIndex,edges_[eIndex][1]),
            labelList(faceHull.size(),-1),
            vertexHull
        )
    );

    // Add this edge to the map.
    map.addEdge(newEdgeIndex);

    // Remove the existing edge from the pointEdges list
    // of the modified point, and add it to the new point
    meshOps::sizeDownList(eIndex, pointEdges_[edges_[eIndex][1]]);
    meshOps::sizeUpList(eIndex, pointEdges_[newPointIndex]);

    // Modify the existing edge
    edges_[eIndex][1] = newPointIndex;

    // Add this edge to the map.
    // Although this edge isn't technically 'added', it's
    // required for coupled patch mapping.
    map.addEdge(eIndex);

    // Keep track of added entities
    labelList addedCellIndices(cellHull.size(),-1);
    labelList addedFaceIndices(faceHull.size(),-1);
    labelList addedEdgeIndices(faceHull.size(),-1);
    labelList addedIntFaceIndices(faceHull.size(),-1);

    // Now loop through the hull and bisect individual entities
    forAll(vertexHull, indexI)
    {
        // Modify the existing face
        meshOps::replaceLabel
        (
            edges_[newEdgeIndex][1],
            newPointIndex,
            faces_[faceHull[indexI]]
        );

        // Modify edgePoints for the edge
        meshOps::replaceLabel
        (
            edges_[newEdgeIndex][1],
            newPointIndex,
            edgePoints_[ringEntities[0][indexI]]
        );

        // Obtain circular indices
        label nextI = vertexHull.fcIndex(indexI);
        label prevI = vertexHull.rcIndex(indexI);

        // Check if this is an interior/boundary face
        if (cellHull[indexI] != -1)
        {
            // Create a new cell. Add it for now, but update later.
            cell newCell(4);

            addedCellIndices[indexI] =
            (
                insertCell(newCell, lengthScale_[cellHull[indexI]])
            );

            // Add this cell to the map.
            map.addCell(addedCellIndices[indexI]);

            // Configure the interior face
            tmpTriFace[0] = vertexHull[nextI];
            tmpTriFace[1] = vertexHull[indexI];
            tmpTriFace[2] = newPointIndex;

            // Insert the face
            addedIntFaceIndices[indexI] =
            (
                insertFace
                (
                    -1,
                    tmpTriFace,
                    cellHull[indexI],
                    addedCellIndices[indexI]
                )
            );

            // Add a faceEdges entry as well
            faceEdges_.append(tmpFaceEdges);

            // Add this face to the map.
            map.addFace(addedIntFaceIndices[indexI]);

            // Add to the new cell
            newCell[0] = addedIntFaceIndices[indexI];

            // Modify the existing ring face connected to newEdge[1]
            label replaceFace = ringEntities[3][indexI];

            // Check if face reversal is necessary
            if (owner_[replaceFace] == cellHull[indexI])
            {
                if (neighbour_[replaceFace] == -1)
                {
                    // Change the owner
                    owner_[replaceFace] = addedCellIndices[indexI];
                }
                else
                {
                    // This face has to be reversed
                    faces_[replaceFace] = faces_[replaceFace].reverseFace();
                    owner_[replaceFace] = neighbour_[replaceFace];
                    neighbour_[replaceFace] = addedCellIndices[indexI];

                    setFlip(replaceFace);
                }
            }
            else
            {
                // Keep owner, but change neighbour
                neighbour_[replaceFace] = addedCellIndices[indexI];
            }

            // Modify the edge on the ring.
            // Add the new interior face to edgeFaces.
            meshOps::sizeUpList
            (
                addedIntFaceIndices[indexI],
                edgeFaces_[edgeHull[indexI]]
            );

            // Insert the new point to edgePoints for the ring edge
            meshOps::insertLabel
            (
                newPointIndex,
                edges_[eIndex][0],
                edges_[newEdgeIndex][1],
                edgePoints_[edgeHull[indexI]]
            );

            // Add this edge to faceEdges for the new interior face
            faceEdges_[addedIntFaceIndices[indexI]][0] = edgeHull[indexI];

            // Replace face labels
            meshOps::replaceLabel
            (
                replaceFace,
                addedIntFaceIndices[indexI],
                cells_[cellHull[indexI]]
            );

            // Add to the new cell
            newCell[1] = replaceFace;

            // Check if this is a boundary face
            if (cellHull[prevI] == -1)
            {
                // Configure the boundary face
                tmpTriFace[0] = newPointIndex;
                tmpTriFace[1] = edges_[newEdgeIndex][1];
                tmpTriFace[2] = vertexHull[indexI];

                // Insert the face
                addedFaceIndices[indexI] =
                (
                    insertFace
                    (
                        whichPatch(faceHull[indexI]),
                        tmpTriFace,
                        addedCellIndices[indexI],
                        -1
                    )
                );

                // Add this face to the map.
                map.addFace(addedFaceIndices[indexI]);

                // Configure edgeFaces
                tmpEdgeFaces[0] = faceHull[indexI];
                tmpEdgeFaces[1] = addedIntFaceIndices[indexI];
                tmpEdgeFaces[2] = addedFaceIndices[indexI];

                // Configure edgePoints
                tmpEdgePoints[0] = edges_[eIndex][0];
                tmpEdgePoints[1] = vertexHull[nextI];
                tmpEdgePoints[2] = edges_[newEdgeIndex][1];

                // Add an edge
                addedEdgeIndices[indexI] =
                (
                    insertEdge
                    (
                        whichPatch(faceHull[indexI]),
                        edge(newPointIndex,vertexHull[indexI]),
                        tmpEdgeFaces,
                        tmpEdgePoints
                    )
                );

                // Add this edge to the map.
                map.addEdge(addedEdgeIndices[indexI]);

                // Add this edge to the interior-face faceEdges entry
                faceEdges_[addedIntFaceIndices[indexI]][1] =
                (
                    addedEdgeIndices[indexI]
                );

                // Configure faceEdges for this boundary face
                tmpFaceEdges[0] = addedEdgeIndices[indexI];
                tmpFaceEdges[1] = newEdgeIndex;
                tmpFaceEdges[2] = ringEntities[2][indexI];

                // Modify faceEdges for the hull face
                meshOps::replaceLabel
                (
                    ringEntities[2][indexI],
                    addedEdgeIndices[indexI],
                    faceEdges_[faceHull[indexI]]
                );

                // Modify edgeFaces for the edge connected to newEdge[1]
                meshOps::replaceLabel
                (
                    faceHull[indexI],
                    addedFaceIndices[indexI],
                    edgeFaces_[ringEntities[2][indexI]]
                );

                // Modify edgePoints for the edge connected to newEdge[1]
                meshOps::replaceLabel
                (
                    edges_[eIndex][0],
                    newPointIndex,
                    edgePoints_[ringEntities[2][indexI]]
                );

                // Add the faceEdges entry
                faceEdges_.append(tmpFaceEdges);

                // Add an entry to edgeFaces
                edgeFaces_[newEdgeIndex][indexI] = addedFaceIndices[indexI];

                // Add an entry for this cell
                newCell[2] = addedFaceIndices[indexI];
            }
            else
            // Check if a cell was added before this
            if (addedCellIndices[prevI] != -1)
            {
                // Configure the interior face
                tmpTriFace[0] = vertexHull[indexI];
                tmpTriFace[1] = edges_[newEdgeIndex][1];
                tmpTriFace[2] = newPointIndex;

                // Insert the face
                addedFaceIndices[indexI] =
                (
                    insertFace
                    (
                        -1,
                        tmpTriFace,
                        addedCellIndices[prevI],
                        addedCellIndices[indexI]
                    )
                );

                // Add this face to the map.
                map.addFace(addedFaceIndices[indexI]);

                // Configure edgeFaces
                tmpIntEdgeFaces[0] = faceHull[indexI];
                tmpIntEdgeFaces[1] = addedIntFaceIndices[indexI];
                tmpIntEdgeFaces[2] = addedFaceIndices[indexI];
                tmpIntEdgeFaces[3] = addedIntFaceIndices[prevI];

                // Configure edgePoints
                tmpIntEdgePoints[0] = edges_[eIndex][0];
                tmpIntEdgePoints[1] = vertexHull[nextI];
                tmpIntEdgePoints[2] = edges_[newEdgeIndex][1];
                tmpIntEdgePoints[3] = vertexHull[prevI];

                // Add an internal edge
                addedEdgeIndices[indexI] =
                (
                    insertEdge
                    (
                        -1,
                        edge(newPointIndex,vertexHull[indexI]),
                        tmpIntEdgeFaces,
                        tmpIntEdgePoints
                    )
                );

                // Add this edge to the map.
                map.addEdge(addedEdgeIndices[indexI]);

                // Add this edge to the interior-face faceEdges entry..
                faceEdges_[addedIntFaceIndices[indexI]][1] =
                (
                    addedEdgeIndices[indexI]
                );

                // ... and to the previous interior face as well
                faceEdges_[addedIntFaceIndices[prevI]][2] =
                (
                    addedEdgeIndices[indexI]
                );

                // Configure faceEdges for this split interior face
                tmpFaceEdges[0] = addedEdgeIndices[indexI];
                tmpFaceEdges[1] = newEdgeIndex;
                tmpFaceEdges[2] = ringEntities[2][indexI];

                // Modify faceEdges for the hull face
                meshOps::replaceLabel
                (
                    ringEntities[2][indexI],
                    addedEdgeIndices[indexI],
                    faceEdges_[faceHull[indexI]]
                );

                // Modify edgeFaces for the edge connected to newEdge[1]
                meshOps::replaceLabel
                (
                    faceHull[indexI],
                    addedFaceIndices[indexI],
                    edgeFaces_[ringEntities[2][indexI]]
                );

                // Modify edgePoints for the edge connected to newEdge[1]
                meshOps::replaceLabel
                (
                    edges_[eIndex][0],
                    newPointIndex,
                    edgePoints_[ringEntities[2][indexI]]
                );

                // Add the faceEdges entry
                faceEdges_.append(tmpFaceEdges);

                // Add an entry to edgeFaces
                edgeFaces_[newEdgeIndex][indexI] = addedFaceIndices[indexI];

                // Add an entry for this cell
                newCell[2] = addedFaceIndices[indexI];

                // Make the final entry for the previous cell
                cells_[addedCellIndices[prevI]][3] = addedFaceIndices[indexI];
            }

            // Do the first interior face at the end
            if (indexI == vertexHull.size() - 1)
            {
                // Configure the interior face
                tmpTriFace[0] = newPointIndex;
                tmpTriFace[1] = edges_[newEdgeIndex][1];
                tmpTriFace[2] = vertexHull[0];

                // Insert the face
                addedFaceIndices[0] =
                (
                    insertFace
                    (
                        -1,
                        tmpTriFace,
                        addedCellIndices[0],
                        addedCellIndices[indexI]
                    )
                );

                // Add this face to the map.
                map.addFace(addedFaceIndices[0]);

                // Configure edgeFaces
                tmpIntEdgeFaces[0] = faceHull[0];
                tmpIntEdgeFaces[1] = addedIntFaceIndices[0];
                tmpIntEdgeFaces[2] = addedFaceIndices[0];
                tmpIntEdgeFaces[3] = addedIntFaceIndices[indexI];

                // Configure edgePoints
                tmpIntEdgePoints[0] = edges_[eIndex][0];
                tmpIntEdgePoints[1] = vertexHull[1];
                tmpIntEdgePoints[2] = edges_[newEdgeIndex][1];
                tmpIntEdgePoints[3] = vertexHull[indexI];

                // Add an internal edge
                addedEdgeIndices[0] =
                (
                    insertEdge
                    (
                        -1,
                        edge(newPointIndex,vertexHull[0]),
                        tmpIntEdgeFaces,
                        tmpIntEdgePoints
                    )
                );

                // Add this edge to the map.
                map.addEdge(addedEdgeIndices[0]);

                // Add this edge to the interior-face faceEdges entry..
                faceEdges_[addedIntFaceIndices[0]][1] =
                (
                    addedEdgeIndices[0]
                );

                // ... and to the previous interior face as well
                faceEdges_[addedIntFaceIndices[indexI]][2] =
                (
                    addedEdgeIndices[0]
                );

                // Configure faceEdges for the first split face
                tmpFaceEdges[0] = addedEdgeIndices[0];
                tmpFaceEdges[1] = newEdgeIndex;
                tmpFaceEdges[2] = ringEntities[2][0];

                // Modify faceEdges for the hull face
                meshOps::replaceLabel
                (
                    ringEntities[2][0],
                    addedEdgeIndices[0],
                    faceEdges_[faceHull[0]]
                );

                // Modify edgeFaces for the edge connected to newEdge[1]
                meshOps::replaceLabel
                (
                    faceHull[0],
                    addedFaceIndices[0],
                    edgeFaces_[ringEntities[2][0]]
                );

                // Modify edgePoints for the edge connected to newEdge[1]
                meshOps::replaceLabel
                (
                    edges_[eIndex][0],
                    newPointIndex,
                    edgePoints_[ringEntities[2][0]]
                );

                // Add the faceEdges entry
                faceEdges_.append(tmpFaceEdges);

                // Add an entry to edgeFaces
                edgeFaces_[newEdgeIndex][0] = addedFaceIndices[0];

                // Add an entry for this cell
                newCell[3] = addedFaceIndices[0];

                // Make the final entry for the first cell
                cells_[addedCellIndices[0]][2] = addedFaceIndices[0];
            }

            // Update the cell list with the new cell.
            cells_[addedCellIndices[indexI]] = newCell;
        }
        else
        {
            // Configure the final boundary face
            tmpTriFace[0] = vertexHull[indexI];
            tmpTriFace[1] = edges_[newEdgeIndex][1];
            tmpTriFace[2] = newPointIndex;

            // Insert the face
            addedFaceIndices[indexI] =
            (
                insertFace
                (
                    whichPatch(faceHull[indexI]),
                    tmpTriFace,
                    addedCellIndices[prevI],
                    -1
                )
            );

            // Add this face to the map.
            map.addFace(addedFaceIndices[indexI]);

            // Configure edgeFaces
            tmpEdgeFaces[0] = addedFaceIndices[indexI];
            tmpEdgeFaces[1] = addedIntFaceIndices[prevI];
            tmpEdgeFaces[2] = faceHull[indexI];

            // Configure edgePoints
            tmpEdgePoints[0] = edges_[newEdgeIndex][1];
            tmpEdgePoints[1] = vertexHull[prevI];
            tmpEdgePoints[2] = edges_[eIndex][0];

            // Add an edge
            addedEdgeIndices[indexI] =
            (
                insertEdge
                (
                    whichPatch(faceHull[indexI]),
                    edge(newPointIndex,vertexHull[indexI]),
                    tmpEdgeFaces,
                    tmpEdgePoints
                )
            );

            // Add this edge to the map.
            map.addEdge(addedEdgeIndices[indexI]);

            // Add a faceEdges entry to the previous interior face
            faceEdges_[addedIntFaceIndices[prevI]][2] =
            (
                addedEdgeIndices[indexI]
            );

            // Configure faceEdges for the final boundary face
            tmpFaceEdges[0] = addedEdgeIndices[indexI];
            tmpFaceEdges[1] = newEdgeIndex;
            tmpFaceEdges[2] = ringEntities[2][indexI];

            // Modify faceEdges for the hull face
            meshOps::replaceLabel
            (
                ringEntities[2][indexI],
                addedEdgeIndices[indexI],
                faceEdges_[faceHull[indexI]]
            );

            // Modify edgeFaces for the edge connected to newEdge[1]
            meshOps::replaceLabel
            (
                faceHull[indexI],
                addedFaceIndices[indexI],
                edgeFaces_[ringEntities[2][indexI]]
            );

            // Modify edgePoints for the edge connected to newEdge[1]
            meshOps::replaceLabel
            (
                edges_[eIndex][0],
                newPointIndex,
                edgePoints_[ringEntities[2][indexI]]
            );

            // Add the faceEdges entry
            faceEdges_.append(tmpFaceEdges);

            // Add an entry to edgeFaces
            edgeFaces_[newEdgeIndex][indexI] = addedFaceIndices[indexI];

            // Make the final entry for the previous cell
            cells_[addedCellIndices[prevI]][3] = addedFaceIndices[indexI];
        }
    }

    // Now that all old / new cells possess correct connectivity,
    // compute mapping information.
    forAll(cellHull, indexI)
    {
        if (cellHull[indexI] == -1)
        {
            continue;
        }

        // Set mapping for both new and modified cells.
        FixedList<label, 2> cmIndex;

        cmIndex[0] = cellHull[indexI];
        cmIndex[1] = addedCellIndices[indexI];

        // Fill-in candidate mapping information
        labelList mC(1, cellHull[indexI]);

        forAll(cmIndex, cmI)
        {
            // Set the mapping for this cell
            setCellMapping(cmIndex[cmI], mC);
        }
    }

    // Set mapping information for old / new faces
    forAll(faceHull, indexI)
    {
        // Interior faces get default mapping
        if (addedIntFaceIndices[indexI] > -1)
        {
            setFaceMapping(addedIntFaceIndices[indexI]);
        }

        // Decide between default / weighted mapping
        // based on boundary information
        if (whichPatch(faceHull[indexI]) == -1)
        {
            // Interior faces get default mapping
            setFaceMapping(faceHull[indexI]);
            setFaceMapping(addedFaceIndices[indexI]);
        }
        else
        {
            // Compute mapping weights for boundary faces
            FixedList<label, 2> fmIndex;

            fmIndex[0] = faceHull[indexI];
            fmIndex[1] = addedFaceIndices[indexI];

            // Fill-in candidate mapping information
            labelList mF(1, faceHull[indexI]);

            forAll(fmIndex, fmI)
            {
                // Set the mapping for this face
                setFaceMapping(fmIndex[fmI], mF);
            }
        }
    }

    if (debug > 2)
    {
        label bPatch = whichEdgePatch(eIndex);

        if (bPatch == -1)
        {
            Info << "Patch: Internal" << endl;
        }
        else
        {
            Info << "Patch: " << boundaryMesh()[bPatch].name() << endl;
        }

        Info << "EdgePoints: " << vertexHull << endl;
        Info << "Edges: " << edgeHull << endl;
        Info << "Faces: " << faceHull << endl;
        Info << "Cells: " << cellHull << endl;

        Info << "Modified cells: " << endl;

        forAll(cellHull, cellI)
        {
            if (cellHull[cellI] == -1)
            {
                continue;
            }

            Info << cellHull[cellI] << ":: "
                 << cells_[cellHull[cellI]]
                 << endl;
        }

        Info << "Added cells: " << endl;

        forAll(addedCellIndices, cellI)
        {
            if (addedCellIndices[cellI] == -1)
            {
                continue;
            }

            Info << addedCellIndices[cellI] << ":: "
                 << cells_[addedCellIndices[cellI]] << nl
                 << "lengthScale: " << lengthScale_[addedCellIndices[cellI]]
                 << endl;
        }

        Info << "Modified faces: " << endl;

        forAll(faceHull, faceI)
        {
            Info << faceHull[faceI] << ":: "
                 << faces_[faceHull[faceI]] << ": "
                 << owner_[faceHull[faceI]] << ": "
                 << neighbour_[faceHull[faceI]] << " "
                 << "faceEdges:: " << faceEdges_[faceHull[faceI]]
                 << endl;
        }

        Info << "Added faces: " << endl;

        forAll(addedFaceIndices, faceI)
        {
            Info << addedFaceIndices[faceI] << ":: "
                 << faces_[addedFaceIndices[faceI]] << ": "
                 << owner_[addedFaceIndices[faceI]] << ": "
                 << neighbour_[addedFaceIndices[faceI]] << " "
                 << "faceEdges:: " << faceEdges_[addedFaceIndices[faceI]]
                 << endl;
        }

        forAll(addedIntFaceIndices, faceI)
        {
            if (addedIntFaceIndices[faceI] == -1)
            {
                continue;
            }

            Info << addedIntFaceIndices[faceI] << ":: "
                 << faces_[addedIntFaceIndices[faceI]] << ": "
                 << owner_[addedIntFaceIndices[faceI]] << ": "
                 << neighbour_[addedIntFaceIndices[faceI]] << " "
                 << "faceEdges:: "
                 << faceEdges_[addedIntFaceIndices[faceI]]
                 << endl;
        }

        Info << "New edge:: " << newEdgeIndex
             << ": " << edges_[newEdgeIndex] << nl
             << " edgeFaces:: " << edgeFaces_[newEdgeIndex] << nl
             << " edgePoints:: " << edgePoints_[newEdgeIndex]
             << endl;

        Info << "Added edges: " << endl;

        forAll(addedEdgeIndices, edgeI)
        {
            Info << addedEdgeIndices[edgeI]
                 << ":: " << edges_[addedEdgeIndices[edgeI]] << nl
                 << " edgeFaces:: " << edgeFaces_[addedEdgeIndices[edgeI]] << nl
                 << " edgePoints:: " << edgePoints_[addedEdgeIndices[edgeI]]
                 << endl;
        }

        Info << "New Point:: " << newPointIndex << endl;
        Info << "pointEdges:: " << pointEdges_[newPointIndex] << endl;

        // Write out VTK files after change
        if (debug > 3)
        {
            labelList newHull(cellHull.size() + addedCellIndices.size(), 0);

            // Combine both lists into one.
            forAll(cellHull, i)
            {
                newHull[i] = cellHull[i];
            }

            label start = cellHull.size();

            for(label i = start; i < newHull.size(); i++)
            {
                newHull[i] = addedCellIndices[i - start];
            }

            writeVTK
            (
                Foam::name(eIndex)
              + "_Bisect_1",
                newHull
            );
        }
    }

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    statistics_[3]++;

    // Increment the number of modifications
    statistics_[0]++;

    // Specify that the operation was successful
    map.type() = 1;

    // Return the changeMap
    return map;
}

// Method for the trisection of a face in 3D
// - Returns a changeMap with a type specifying:
//     1: Trisection was successful
//    -1: Trisection failed since max number of topo-changes was reached.
//    -2: Trisection failed since resulting quality would be really bad.
// - AddedPoint is the index of the newly added point.
const changeMap dynamicTopoFvMesh::trisectFace
(
    const label fIndex,
    bool checkOnly,
    bool forceOp
)
{
    // Face trisection performs the following operations:
    //      [1] Add a point at middle of the face
    //      [2] Remove the face and add three new faces in place.
    //      [3] Add three cells for each trisected cell (remove the originals).
    //      [4] Create one internal edge for each trisected cell.
    //      [5] Create three edges for the trisected face.
    //      [6] Create three internal faces for each trisected cell.
    //      Update faceEdges, edgeFaces and edgePoints information.

    // Figure out which thread this is...
    label tIndex = self();

    // Prepare the changeMaps
    changeMap map, slaveMap;

    if
    (
        (statistics_[0] > maxModifications_) &&
        (maxModifications_ > -1)
    )
    {
        // Reached the max allowable topo-changes.
        stack(tIndex).clear();

        return map;
    }

    // Sanity check: Is the index legitimate?
    if (fIndex < 0)
    {
        FatalErrorIn
        (
            "const changeMap dynamicTopoFvMesh::trisectFace\n"
            "(\n"
            "    const label fIndex,\n"
            "    bool checkOnly,\n"
            "    bool forceOp\n"
            ")\n"
        )
            << " Invalid index: " << fIndex
            << abort(FatalError);
    }

    // Before we trisect this face, check whether the operation will
    // yield an acceptable cell-quality.
    scalar minQ = 0.0;

    if ((minQ = computeTrisectionQuality(fIndex)) < sliverThreshold_)
    {
        // Check if the quality is actually valid before forcing it.
        if (forceOp && (minQ < 0.0))
        {
            FatalErrorIn
            (
                "const changeMap dynamicTopoFvMesh::trisectFace\n"
                "(\n"
                "    const label fIndex,\n"
                "    bool checkOnly,\n"
                "    bool forceOp\n"
                ")\n"
            )
                << " Forcing trisection on face: " << fIndex
                << " will yield an invalid cell."
                << abort(FatalError);
        }
        else
        if (!forceOp)
        {
            map.type() = -2;
            return map;
        }
    }

    // Are we performing only checks?
    if (checkOnly)
    {
        map.type() = 1;
        return map;
    }

    // Update number of surface bisections, if necessary.
    if (whichPatch(fIndex) > -1)
    {
        statistics_[5]++;
    }

    // Hull variables
    face tmpTriFace(3);
    labelList newTriEdgeFaces(3), newTriEdgePoints(3);
    labelList newQuadEdgeFaces(4), newQuadEdgePoints(4);

    FixedList<label,2> apexPoint(-1);
    FixedList<face, 3> checkFace(face(3));
    FixedList<label,5> newEdgeIndex(-1);
    FixedList<label,9> newFaceIndex(-1);
    FixedList<label,6> newCellIndex(-1);
    FixedList<cell, 6> newTetCell(cell(4));
    FixedList<labelList, 9> newFaceEdges(labelList(3));

    // Counters for entities
    FixedList<label, 9> nE(0);
    FixedList<label, 6> nF(0);

    // Determine the two cells to be removed
    FixedList<label,2> cellsForRemoval;
    cellsForRemoval[0] = owner_[fIndex];
    cellsForRemoval[1] = neighbour_[fIndex];

    if (debug > 1)
    {
        Info << nl << nl
             << "Face: " << fIndex
             << ": " << faces_[fIndex]
             << " is to be trisected. " << endl;

        // Write out VTK files prior to change
        if (debug > 3)
        {
            labelList vtkCells;

            if (neighbour_[fIndex] == -1)
            {
                vtkCells.setSize(1);
                vtkCells[0] = owner_[fIndex];
            }
            else
            {
                vtkCells.setSize(2);
                vtkCells[0] = owner_[fIndex];
                vtkCells[1] = neighbour_[fIndex];
            }

            writeVTK
            (
                Foam::name(fIndex)
              + "Trisect_0",
                vtkCells
            );
        }
    }

    labelList mP(3, -1);

    // Fill in mapping information
    mP[0] = faces_[fIndex][0];
    mP[1] = faces_[fIndex][1];
    mP[2] = faces_[fIndex][2];

    // Add a new point to the end of the list
    scalar oT = (1.0/3.0);

    label newPointIndex =
    (
        insertPoint
        (
            oT * (points_[mP[0]] + points_[mP[1]] + points_[mP[2]]),
            oT * (oldPoints_[mP[0]] + oldPoints_[mP[1]] + oldPoints_[mP[2]]),
            mP
        )
    );

    // Add this point to the map.
    map.addPoint(newPointIndex);

    // Add three new cells to the end of the cell list
    for (label i = 0; i < 3; i++)
    {
        scalar parentScale = -1.0;

        if (edgeRefinement_)
        {
            parentScale = lengthScale_[cellsForRemoval[0]];
        }

        newCellIndex[i] = insertCell(newTetCell[i], parentScale);

        // Add cells to the map
        map.addCell(newCellIndex[i]);
    }

    // Find the apex point for this cell
    apexPoint[0] =
    (
        meshOps::tetApexPoint
        (
            owner_[fIndex],
            fIndex,
            faces_,
            cells_
        )
    );

    // Insert three new internal faces

    // First face: Owner: newCellIndex[0], Neighbour: newCellIndex[1]
    tmpTriFace[0] = newPointIndex;
    tmpTriFace[1] = faces_[fIndex][0];
    tmpTriFace[2] = apexPoint[0];

    newFaceIndex[0] =
    (
        insertFace
        (
            -1,
            tmpTriFace,
            newCellIndex[0],
            newCellIndex[1]
        )
    );

    // Second face: Owner: newCellIndex[1], Neighbour: newCellIndex[2]
    tmpTriFace[0] = newPointIndex;
    tmpTriFace[1] = faces_[fIndex][1];
    tmpTriFace[2] = apexPoint[0];

    newFaceIndex[1] =
    (
        insertFace
        (
            -1,
            tmpTriFace,
            newCellIndex[1],
            newCellIndex[2]
        )
    );

    // Third face: Owner: newCellIndex[0], Neighbour: newCellIndex[2]
    tmpTriFace[0] = newPointIndex;
    tmpTriFace[1] = apexPoint[0];
    tmpTriFace[2] = faces_[fIndex][2];

    newFaceIndex[2] =
    (
        insertFace
        (
            -1,
            tmpTriFace,
            newCellIndex[0],
            newCellIndex[2]
        )
    );

    // Add an entry to edgeFaces
    newTriEdgeFaces[0] = newFaceIndex[0];
    newTriEdgeFaces[1] = newFaceIndex[1];
    newTriEdgeFaces[2] = newFaceIndex[2];

    // Add an entry for edgePoints as well
    newTriEdgePoints[0] = faces_[fIndex][0];
    newTriEdgePoints[1] = faces_[fIndex][1];
    newTriEdgePoints[2] = faces_[fIndex][2];

    // Add a new internal edge to the mesh
    newEdgeIndex[0] =
    (
        insertEdge
        (
            -1,
            edge
            (
               newPointIndex,
               apexPoint[0]
            ),
            newTriEdgeFaces,
            newTriEdgePoints
        )
    );

    // Configure faceEdges with the new internal edge
    newFaceEdges[0][nE[0]++] = newEdgeIndex[0];
    newFaceEdges[1][nE[1]++] = newEdgeIndex[0];
    newFaceEdges[2][nE[2]++] = newEdgeIndex[0];

    // Add the newly created faces to cells
    newTetCell[0][nF[0]++] = newFaceIndex[0];
    newTetCell[0][nF[0]++] = newFaceIndex[2];
    newTetCell[1][nF[1]++] = newFaceIndex[0];
    newTetCell[1][nF[1]++] = newFaceIndex[1];
    newTetCell[2][nF[2]++] = newFaceIndex[1];
    newTetCell[2][nF[2]++] = newFaceIndex[2];

    // Define the three faces to check for orientation:
    checkFace[0][0] = faces_[fIndex][2];
    checkFace[0][1] = apexPoint[0];
    checkFace[0][2] = faces_[fIndex][0];

    checkFace[1][0] = faces_[fIndex][0];
    checkFace[1][1] = apexPoint[0];
    checkFace[1][2] = faces_[fIndex][1];

    checkFace[2][0] = faces_[fIndex][1];
    checkFace[2][1] = apexPoint[0];
    checkFace[2][2] = faces_[fIndex][2];

    // Check the orientation of faces on the first cell.
    forAll(cells_[owner_[fIndex]], faceI)
    {
        label faceIndex = cells_[owner_[fIndex]][faceI];

        if (faceIndex == fIndex)
        {
            continue;
        }

        const face& faceToCheck = faces_[faceIndex];
        label cellIndex = cellsForRemoval[0];
        label newIndex = -1;

        // Check against faces.
        if (triFace::compare(triFace(faceToCheck), triFace(checkFace[0])))
        {
            newIndex = newCellIndex[0];
            newTetCell[0][nF[0]++] = faceIndex;
        }
        else
        if (triFace::compare(triFace(faceToCheck), triFace(checkFace[1])))
        {
            newIndex = newCellIndex[1];
            newTetCell[1][nF[1]++] = faceIndex;
        }
        else
        if (triFace::compare(triFace(faceToCheck), triFace(checkFace[2])))
        {
            newIndex = newCellIndex[2];
            newTetCell[2][nF[2]++] = faceIndex;
        }
        else
        {
            // Something's terribly wrong.
            FatalErrorIn
            (
                "const changeMap dynamicTopoFvMesh::trisectFace\n"
                "(\n"
                "    const label fIndex,\n"
                "    bool checkOnly,\n"
                "    bool forceOp\n"
                ")\n"
            )
                << "Failed to determine a face match."
                << abort(FatalError);
        }

        // Check if a face-flip is necessary
        if (owner_[faceIndex] == cellIndex)
        {
            if (neighbour_[faceIndex] == -1)
            {
                // Change the owner
                owner_[faceIndex] = newIndex;
            }
            else
            {
                // Flip this face
                faces_[faceIndex] = faceToCheck.reverseFace();
                owner_[faceIndex] = neighbour_[faceIndex];
                neighbour_[faceIndex] = newIndex;

                setFlip(faceIndex);
            }
        }
        else
        {
            // Flip is unnecessary. Just update neighbour
            neighbour_[faceIndex] = newIndex;
        }
    }

    if (cellsForRemoval[1] == -1)
    {
        // Boundary face. Determine its patch.
        label facePatch = whichPatch(fIndex);

        // Add three new boundary faces.

        // Fourth face: Owner: newCellIndex[0], Neighbour: -1
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = faces_[fIndex][2];
        tmpTriFace[2] = faces_[fIndex][0];

        newFaceIndex[3] =
        (
            insertFace
            (
                facePatch,
                tmpTriFace,
                newCellIndex[0],
                -1
            )
        );

        // Fifth face: Owner: newCellIndex[1], Neighbour: -1
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = faces_[fIndex][0];
        tmpTriFace[2] = faces_[fIndex][1];

        newFaceIndex[4] =
        (
            insertFace
            (
                facePatch,
                tmpTriFace,
                newCellIndex[1],
                -1
            )
        );

        // Sixth face: Owner: newCellIndex[2], Neighbour: -1
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = faces_[fIndex][1];
        tmpTriFace[2] = faces_[fIndex][2];

        newFaceIndex[5] =
        (
            insertFace
            (
                facePatch,
                tmpTriFace,
                newCellIndex[2],
                -1
            )
        );

        // Add the newly created faces to cells
        newTetCell[0][nF[0]++] = newFaceIndex[3];
        newTetCell[1][nF[1]++] = newFaceIndex[4];
        newTetCell[2][nF[2]++] = newFaceIndex[5];

        // Configure edgeFaces and edgePoints for three new boundary edges.
        newTriEdgeFaces[0] = newFaceIndex[4];
        newTriEdgeFaces[1] = newFaceIndex[0];
        newTriEdgeFaces[2] = newFaceIndex[3];

        newTriEdgePoints[0] = faces_[fIndex][1];
        newTriEdgePoints[1] = apexPoint[0];
        newTriEdgePoints[2] = faces_[fIndex][2];

        newEdgeIndex[1] =
        (
            insertEdge
            (
                facePatch,
                edge
                (
                   newPointIndex,
                   faces_[fIndex][0]
                ),
                newTriEdgeFaces,
                newTriEdgePoints
            )
        );

        newTriEdgeFaces[0] = newFaceIndex[5];
        newTriEdgeFaces[1] = newFaceIndex[1];
        newTriEdgeFaces[2] = newFaceIndex[4];

        newTriEdgePoints[0] = faces_[fIndex][2];
        newTriEdgePoints[1] = apexPoint[0];
        newTriEdgePoints[2] = faces_[fIndex][0];

        newEdgeIndex[2] =
        (
            insertEdge
            (
                facePatch,
                edge
                (
                   newPointIndex,
                   faces_[fIndex][1]
                ),
                newTriEdgeFaces,
                newTriEdgePoints
            )
        );

        newTriEdgeFaces[0] = newFaceIndex[3];
        newTriEdgeFaces[1] = newFaceIndex[2];
        newTriEdgeFaces[2] = newFaceIndex[5];

        newTriEdgePoints[0] = faces_[fIndex][0];
        newTriEdgePoints[1] = apexPoint[0];
        newTriEdgePoints[2] = faces_[fIndex][1];

        newEdgeIndex[3] =
        (
            insertEdge
            (
                facePatch,
                edge
                (
                   newPointIndex,
                   faces_[fIndex][2]
                ),
                newTriEdgeFaces,
                newTriEdgePoints
            )
        );

        // Configure faceEdges with the three new edges.
        newFaceEdges[0][nE[0]++] = newEdgeIndex[1];
        newFaceEdges[1][nE[1]++] = newEdgeIndex[2];
        newFaceEdges[2][nE[2]++] = newEdgeIndex[3];

        newFaceEdges[3][nE[3]++] = newEdgeIndex[1];
        newFaceEdges[3][nE[3]++] = newEdgeIndex[3];
        newFaceEdges[4][nE[4]++] = newEdgeIndex[1];
        newFaceEdges[4][nE[4]++] = newEdgeIndex[2];
        newFaceEdges[5][nE[5]++] = newEdgeIndex[2];
        newFaceEdges[5][nE[5]++] = newEdgeIndex[3];

        // Define the six edges to check while building faceEdges:
        FixedList<edge,6> check;

        check[0][0] = apexPoint[0]; check[0][1] = faces_[fIndex][0];
        check[1][0] = apexPoint[0]; check[1][1] = faces_[fIndex][1];
        check[2][0] = apexPoint[0]; check[2][1] = faces_[fIndex][2];

        check[3][0] = faces_[fIndex][2]; check[3][1] = faces_[fIndex][0];
        check[4][0] = faces_[fIndex][0]; check[4][1] = faces_[fIndex][1];
        check[5][0] = faces_[fIndex][1]; check[5][1] = faces_[fIndex][2];

        // Build a list of cellEdges
        labelHashSet cellEdges;

        forAll(cells_[owner_[fIndex]], faceI)
        {
            const labelList& fEdges =
            (
                faceEdges_[cells_[owner_[fIndex]][faceI]]
            );

            forAll(fEdges, edgeI)
            {
                if (!cellEdges.found(fEdges[edgeI]))
                {
                    cellEdges.insert(fEdges[edgeI]);
                }
            }
        }

        // Loop through cellEdges, and perform appropriate actions.
        forAllIter(labelHashSet, cellEdges, eIter)
        {
            const edge& edgeToCheck = edges_[eIter.key()];

            // Check against the specified edges.
            if (edgeToCheck == check[0])
            {
                meshOps::insertLabel
                (
                    newPointIndex,
                    faces_[fIndex][1],
                    faces_[fIndex][2],
                    edgePoints_[eIter.key()]
                );

                meshOps::sizeUpList(newFaceIndex[0], edgeFaces_[eIter.key()]);
                newFaceEdges[0][nE[0]++] = eIter.key();
            }

            if (edgeToCheck == check[1])
            {
                meshOps::insertLabel
                (
                    newPointIndex,
                    faces_[fIndex][0],
                    faces_[fIndex][2],
                    edgePoints_[eIter.key()]
                );

                meshOps::sizeUpList(newFaceIndex[1], edgeFaces_[eIter.key()]);
                newFaceEdges[1][nE[1]++] = eIter.key();
            }

            if (edgeToCheck == check[2])
            {
                meshOps::insertLabel
                (
                    newPointIndex,
                    faces_[fIndex][0],
                    faces_[fIndex][1],
                    edgePoints_[eIter.key()]
                );

                meshOps::sizeUpList(newFaceIndex[2], edgeFaces_[eIter.key()]);
                newFaceEdges[2][nE[2]++] = eIter.key();
            }

            if (edgeToCheck == check[3])
            {
                meshOps::replaceLabel
                (
                    faces_[fIndex][1],
                    newPointIndex,
                    edgePoints_[eIter.key()]
                );

                meshOps::replaceLabel
                (
                    fIndex,
                    newFaceIndex[3],
                    edgeFaces_[eIter.key()]
                );

                newFaceEdges[3][nE[3]++] = eIter.key();
            }

            if (edgeToCheck == check[4])
            {
                meshOps::replaceLabel
                (
                    faces_[fIndex][2],
                    newPointIndex,
                    edgePoints_[eIter.key()]
                );

                meshOps::replaceLabel
                (
                    fIndex,
                    newFaceIndex[4],
                    edgeFaces_[eIter.key()]
                );

                newFaceEdges[4][nE[4]++] = eIter.key();
            }

            if (edgeToCheck == check[5])
            {
                meshOps::replaceLabel
                (
                    faces_[fIndex][0],
                    newPointIndex,
                    edgePoints_[eIter.key()]
                );

                meshOps::replaceLabel
                (
                    fIndex,
                    newFaceIndex[5],
                    edgeFaces_[eIter.key()]
                );

                newFaceEdges[5][nE[5]++] = eIter.key();
            }
        }

        // Now that faceEdges has been configured, append them to the list.
        for (label i = 0; i < 6; i++)
        {
            faceEdges_.append(newFaceEdges[i]);

            // Add faces to the map.
            map.addFace(newFaceIndex[i]);
        }
    }
    else
    {
        // Add three new cells to the end of the cell list
        for (label i = 3; i < 6; i++)
        {
            scalar parentScale = -1.0;

            if (edgeRefinement_)
            {
                parentScale = lengthScale_[cellsForRemoval[1]];
            }

            newCellIndex[i] = insertCell(newTetCell[i], parentScale);

            // Add to the map.
            map.addCell(newCellIndex[i]);
        }

        // Find the apex point for this cell
        apexPoint[1] =
        (
            meshOps::tetApexPoint
            (
                neighbour_[fIndex],
                fIndex,
                faces_,
                cells_
            )
        );

        // Add six new interior faces.

        // Fourth face: Owner: newCellIndex[0], Neighbour: newCellIndex[3]
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = faces_[fIndex][2];
        tmpTriFace[2] = faces_[fIndex][0];

        newFaceIndex[3] =
        (
            insertFace
            (
                -1,
                tmpTriFace,
                newCellIndex[0],
                newCellIndex[3]
            )
        );

        // Fifth face: Owner: newCellIndex[1], Neighbour: newCellIndex[4]
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = faces_[fIndex][0];
        tmpTriFace[2] = faces_[fIndex][1];

        newFaceIndex[4] =
        (
            insertFace
            (
                -1,
                tmpTriFace,
                newCellIndex[1],
                newCellIndex[4]
            )
        );

        // Sixth face: Owner: newCellIndex[2], Neighbour: newCellIndex[5]
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = faces_[fIndex][1];
        tmpTriFace[2] = faces_[fIndex][2];

        newFaceIndex[5] =
        (
            insertFace
            (
                -1,
                tmpTriFace,
                newCellIndex[2],
                newCellIndex[5]
            )
        );

        // Seventh face: Owner: newCellIndex[3], Neighbour: newCellIndex[4]
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = apexPoint[1];
        tmpTriFace[2] = faces_[fIndex][0];

        newFaceIndex[6] =
        (
            insertFace
            (
                -1,
                tmpTriFace,
                newCellIndex[3],
                newCellIndex[4]
            )
        );

        // Eighth face: Owner: newCellIndex[4], Neighbour: newCellIndex[5]
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = apexPoint[1];
        tmpTriFace[2] = faces_[fIndex][1];

        newFaceIndex[7] =
        (
            insertFace
            (
                -1,
                tmpTriFace,
                newCellIndex[4],
                newCellIndex[5]
            )
        );

        // Ninth face: Owner: newCellIndex[3], Neighbour: newCellIndex[5]
        tmpTriFace[0] = newPointIndex;
        tmpTriFace[1] = faces_[fIndex][2];
        tmpTriFace[2] = apexPoint[1];

        newFaceIndex[8] =
        (
            insertFace
            (
                -1,
                tmpTriFace,
                newCellIndex[3],
                newCellIndex[5]
            )
        );

        // Add the newly created faces to cells
        newTetCell[3][nF[3]++] = newFaceIndex[6];
        newTetCell[3][nF[3]++] = newFaceIndex[8];
        newTetCell[4][nF[4]++] = newFaceIndex[6];
        newTetCell[4][nF[4]++] = newFaceIndex[7];
        newTetCell[5][nF[5]++] = newFaceIndex[7];
        newTetCell[5][nF[5]++] = newFaceIndex[8];

        newTetCell[0][nF[0]++] = newFaceIndex[3];
        newTetCell[1][nF[1]++] = newFaceIndex[4];
        newTetCell[2][nF[2]++] = newFaceIndex[5];

        newTetCell[3][nF[3]++] = newFaceIndex[3];
        newTetCell[4][nF[4]++] = newFaceIndex[4];
        newTetCell[5][nF[5]++] = newFaceIndex[5];

        // Define the three faces to check for orientation:
        checkFace[0][0] = faces_[fIndex][2];
        checkFace[0][1] = apexPoint[1];
        checkFace[0][2] = faces_[fIndex][0];

        checkFace[1][0] = faces_[fIndex][0];
        checkFace[1][1] = apexPoint[1];
        checkFace[1][2] = faces_[fIndex][1];

        checkFace[2][0] = faces_[fIndex][1];
        checkFace[2][1] = apexPoint[1];
        checkFace[2][2] = faces_[fIndex][2];

        // Check the orientation of faces on the second cell.
        forAll(cells_[neighbour_[fIndex]], faceI)
        {
            label faceIndex = cells_[neighbour_[fIndex]][faceI];

            if (faceIndex == fIndex)
            {
                continue;
            }

            const face& faceToCheck = faces_[faceIndex];
            label cellIndex = cellsForRemoval[1];
            label newIndex = -1;

            // Check against faces.
            if (triFace::compare(triFace(faceToCheck), triFace(checkFace[0])))
            {
                newIndex = newCellIndex[3];
                newTetCell[3][nF[3]++] = faceIndex;
            }
            else
            if (triFace::compare(triFace(faceToCheck), triFace(checkFace[1])))
            {
                newIndex = newCellIndex[4];
                newTetCell[4][nF[4]++] = faceIndex;
            }
            else
            if (triFace::compare(triFace(faceToCheck), triFace(checkFace[2])))
            {
                newIndex = newCellIndex[5];
                newTetCell[5][nF[5]++] = faceIndex;
            }
            else
            {
                // Something's terribly wrong.
                FatalErrorIn
                (
                    "const changeMap dynamicTopoFvMesh::trisectFace\n"
                    "(\n"
                    "    const label fIndex,\n"
                    "    bool checkOnly,\n"
                    "    bool forceOp\n"
                    ")\n"
                )
                    << "Failed to determine a face match."
                    << abort(FatalError);
            }

            // Check if a face-flip is necessary
            if (owner_[faceIndex] == cellIndex)
            {
                if (neighbour_[faceIndex] == -1)
                {
                    // Change the owner
                    owner_[faceIndex] = newIndex;
                }
                else
                {
                    // Flip this face
                    faces_[faceIndex] = faceToCheck.reverseFace();
                    owner_[faceIndex] = neighbour_[faceIndex];
                    neighbour_[faceIndex] = newIndex;

                    setFlip(faceIndex);
                }
            }
            else
            {
                // Flip is unnecessary. Just update neighbour
                neighbour_[faceIndex] = newIndex;
            }
        }

        // Configure edgeFaces and edgePoints for four new interior edges.
        newQuadEdgeFaces[0] = newFaceIndex[4];
        newQuadEdgeFaces[1] = newFaceIndex[0];
        newQuadEdgeFaces[2] = newFaceIndex[3];
        newQuadEdgeFaces[3] = newFaceIndex[6];

        newQuadEdgePoints[0] = faces_[fIndex][1];
        newQuadEdgePoints[1] = apexPoint[0];
        newQuadEdgePoints[2] = faces_[fIndex][2];
        newQuadEdgePoints[3] = apexPoint[1];

        newEdgeIndex[1] =
        (
            insertEdge
            (
                -1,
                edge
                (
                   newPointIndex,
                   faces_[fIndex][0]
                ),
                newQuadEdgeFaces,
                newQuadEdgePoints
            )
        );

        newQuadEdgeFaces[0] = newFaceIndex[5];
        newQuadEdgeFaces[1] = newFaceIndex[1];
        newQuadEdgeFaces[2] = newFaceIndex[4];
        newQuadEdgeFaces[3] = newFaceIndex[7];

        newQuadEdgePoints[0] = faces_[fIndex][2];
        newQuadEdgePoints[1] = apexPoint[0];
        newQuadEdgePoints[2] = faces_[fIndex][0];
        newQuadEdgePoints[3] = apexPoint[1];

        newEdgeIndex[2] =
        (
            insertEdge
            (
                -1,
                edge
                (
                   newPointIndex,
                   faces_[fIndex][1]
                ),
                newQuadEdgeFaces,
                newQuadEdgePoints
            )
        );

        newQuadEdgeFaces[0] = newFaceIndex[3];
        newQuadEdgeFaces[1] = newFaceIndex[2];
        newQuadEdgeFaces[2] = newFaceIndex[5];
        newQuadEdgeFaces[3] = newFaceIndex[8];

        newQuadEdgePoints[0] = faces_[fIndex][0];
        newQuadEdgePoints[1] = apexPoint[0];
        newQuadEdgePoints[2] = faces_[fIndex][1];
        newQuadEdgePoints[3] = apexPoint[1];

        newEdgeIndex[3] =
        (
            insertEdge
            (
                -1,
                edge
                (
                   newPointIndex,
                   faces_[fIndex][2]
                ),
                newQuadEdgeFaces,
                newQuadEdgePoints
            )
        );

        newTriEdgeFaces[0] = newFaceIndex[6];
        newTriEdgeFaces[1] = newFaceIndex[7];
        newTriEdgeFaces[2] = newFaceIndex[8];

        newTriEdgePoints[0] = faces_[fIndex][0];
        newTriEdgePoints[1] = faces_[fIndex][1];
        newTriEdgePoints[2] = faces_[fIndex][2];

        newEdgeIndex[4] =
        (
            insertEdge
            (
                -1,
                edge
                (
                   apexPoint[1],
                   newPointIndex
                ),
                newTriEdgeFaces,
                newTriEdgePoints
            )
        );

        // Configure faceEdges with the new internal edges
        newFaceEdges[0][nE[0]++] = newEdgeIndex[1];
        newFaceEdges[1][nE[1]++] = newEdgeIndex[2];
        newFaceEdges[2][nE[2]++] = newEdgeIndex[3];

        newFaceEdges[3][nE[3]++] = newEdgeIndex[1];
        newFaceEdges[3][nE[3]++] = newEdgeIndex[3];
        newFaceEdges[4][nE[4]++] = newEdgeIndex[1];
        newFaceEdges[4][nE[4]++] = newEdgeIndex[2];
        newFaceEdges[5][nE[5]++] = newEdgeIndex[2];
        newFaceEdges[5][nE[5]++] = newEdgeIndex[3];

        newFaceEdges[6][nE[6]++] = newEdgeIndex[1];
        newFaceEdges[7][nE[7]++] = newEdgeIndex[2];
        newFaceEdges[8][nE[8]++] = newEdgeIndex[3];

        newFaceEdges[6][nE[6]++] = newEdgeIndex[4];
        newFaceEdges[7][nE[7]++] = newEdgeIndex[4];
        newFaceEdges[8][nE[8]++] = newEdgeIndex[4];

        // Define the nine edges to check while building faceEdges:
        FixedList<edge,9> check;

        check[0][0] = apexPoint[0]; check[0][1] = faces_[fIndex][0];
        check[1][0] = apexPoint[0]; check[1][1] = faces_[fIndex][1];
        check[2][0] = apexPoint[0]; check[2][1] = faces_[fIndex][2];

        check[3][0] = faces_[fIndex][2]; check[3][1] = faces_[fIndex][0];
        check[4][0] = faces_[fIndex][0]; check[4][1] = faces_[fIndex][1];
        check[5][0] = faces_[fIndex][1]; check[5][1] = faces_[fIndex][2];

        check[6][0] = apexPoint[1]; check[6][1] = faces_[fIndex][0];
        check[7][0] = apexPoint[1]; check[7][1] = faces_[fIndex][1];
        check[8][0] = apexPoint[1]; check[8][1] = faces_[fIndex][2];

        // Build a list of cellEdges
        labelHashSet cellEdges;

        forAll(cells_[owner_[fIndex]], faceI)
        {
            const labelList& fEdges =
            (
                faceEdges_[cells_[owner_[fIndex]][faceI]]
            );

            forAll(fEdges, edgeI)
            {
                if (!cellEdges.found(fEdges[edgeI]))
                {
                    cellEdges.insert(fEdges[edgeI]);
                }
            }
        }

        forAll(cells_[neighbour_[fIndex]], faceI)
        {
            const labelList& fEdges =
            (
                faceEdges_[cells_[neighbour_[fIndex]][faceI]]
            );

            forAll(fEdges, edgeI)
            {
                if (!cellEdges.found(fEdges[edgeI]))
                {
                    cellEdges.insert(fEdges[edgeI]);
                }
            }
        }

        // Loop through cellEdges, and perform appropriate actions.
        forAllIter(labelHashSet, cellEdges, eIter)
        {
            const edge& edgeToCheck = edges_[eIter.key()];

            // Check against the specified edges.
            if (edgeToCheck == check[0])
            {
                meshOps::insertLabel
                (
                    newPointIndex,
                    faces_[fIndex][1],
                    faces_[fIndex][2],
                    edgePoints_[eIter.key()]
                );

                meshOps::sizeUpList(newFaceIndex[0], edgeFaces_[eIter.key()]);
                newFaceEdges[0][nE[0]++] = eIter.key();
            }

            if (edgeToCheck == check[1])
            {
                meshOps::insertLabel
                (
                    newPointIndex,
                    faces_[fIndex][0],
                    faces_[fIndex][2],
                    edgePoints_[eIter.key()]
                );

                meshOps::sizeUpList(newFaceIndex[1], edgeFaces_[eIter.key()]);
                newFaceEdges[1][nE[1]++] = eIter.key();
            }

            if (edgeToCheck == check[2])
            {
                meshOps::insertLabel
                (
                    newPointIndex,
                    faces_[fIndex][0],
                    faces_[fIndex][1],
                    edgePoints_[eIter.key()]
                );

                meshOps::sizeUpList(newFaceIndex[2], edgeFaces_[eIter.key()]);
                newFaceEdges[2][nE[2]++] = eIter.key();
            }

            if (edgeToCheck == check[3])
            {
                meshOps::replaceLabel
                (
                    faces_[fIndex][1],
                    newPointIndex,
                    edgePoints_[eIter.key()]
                );

                meshOps::replaceLabel
                (
                    fIndex,
                    newFaceIndex[3],
                    edgeFaces_[eIter.key()]
                );

                newFaceEdges[3][nE[3]++] = eIter.key();
            }

            if (edgeToCheck == check[4])
            {
                meshOps::replaceLabel
                (
                    faces_[fIndex][2],
                    newPointIndex,
                    edgePoints_[eIter.key()]
                );

                meshOps::replaceLabel
                (
                    fIndex,
                    newFaceIndex[4],
                    edgeFaces_[eIter.key()]
                );

                newFaceEdges[4][nE[4]++] = eIter.key();
            }

            if (edgeToCheck == check[5])
            {
                meshOps::replaceLabel
                (
                    faces_[fIndex][0],
                    newPointIndex,
                    edgePoints_[eIter.key()]
                );

                meshOps::replaceLabel
                (
                    fIndex,
                    newFaceIndex[5],
                    edgeFaces_[eIter.key()]
                );

                newFaceEdges[5][nE[5]++] = eIter.key();
            }

            if (edgeToCheck == check[6])
            {
                meshOps::insertLabel
                (
                    newPointIndex,
                    faces_[fIndex][1],
                    faces_[fIndex][2],
                    edgePoints_[eIter.key()]
                );

                meshOps::sizeUpList(newFaceIndex[6], edgeFaces_[eIter.key()]);
                newFaceEdges[6][nE[6]++] = eIter.key();
            }

            if (edgeToCheck == check[7])
            {
                meshOps::insertLabel
                (
                    newPointIndex,
                    faces_[fIndex][0],
                    faces_[fIndex][2],
                    edgePoints_[eIter.key()]
                );

                meshOps::sizeUpList(newFaceIndex[7], edgeFaces_[eIter.key()]);
                newFaceEdges[7][nE[7]++] = eIter.key();
            }

            if (edgeToCheck == check[8])
            {
                meshOps::insertLabel
                (
                    newPointIndex,
                    faces_[fIndex][0],
                    faces_[fIndex][1],
                    edgePoints_[eIter.key()]
                );

                meshOps::sizeUpList(newFaceIndex[8], edgeFaces_[eIter.key()]);
                newFaceEdges[8][nE[8]++] = eIter.key();
            }
        }

        // Now that faceEdges has been configured, append them to the list.
        for (label i = 0; i < 9; i++)
        {
            faceEdges_.append(newFaceEdges[i]);

            // Add faces to the map.
            map.addFace(newFaceIndex[i]);
        }
    }

    // Added edges are those connected to the new point
    const labelList& pointEdges = pointEdges_[newPointIndex];

    forAll(pointEdges, edgeI)
    {
        map.addEdge(pointEdges[edgeI]);
    }

    // Now generate mapping info and remove entities.
    forAll(cellsForRemoval, cellI)
    {
        label cIndex = cellsForRemoval[cellI];

        if (cIndex == -1)
        {
            continue;
        }

        // Fill-in mapping information
        labelList mC(1, cellsForRemoval[cellI]);

        if (cellI == 0)
        {
            for (label i = 0; i < 3; i++)
            {
                // Update the cell list with newly configured cells.
                cells_[newCellIndex[i]] = newTetCell[i];

                setCellMapping(newCellIndex[i], mC);
            }
        }
        else
        {
            for (label i = 3; i < 6; i++)
            {
                // Update the cell list with newly configured cells.
                cells_[newCellIndex[i]] = newTetCell[i];

                setCellMapping(newCellIndex[i], mC);
            }
        }

        removeCell(cIndex);
    }

    // Set default mapping for interior faces.
    for (label i = 0; i < 3; i++)
    {
        setFaceMapping(newFaceIndex[i]);
    }

    if (cellsForRemoval[1] == -1)
    {
        // Set mapping for boundary faces.
        for (label i = 3; i < 6; i++)
        {
            setFaceMapping(newFaceIndex[i], labelList(1, fIndex));
        }
    }
    else
    {
        // Set default mapping for interior faces.
        for (label i = 3; i < 9; i++)
        {
            setFaceMapping(newFaceIndex[i]);
        }
    }

    // Now finally remove the face...
    removeFace(fIndex);

    if (debug > 2)
    {
        Info << "New Point:: " << newPointIndex << endl;

        const labelList& pEdges = pointEdges_[newPointIndex];

        Info << "pointEdges:: " << pEdges << endl;

        Info << "Added edges: " << endl;
        forAll(pEdges, edgeI)
        {
            Info << pEdges[edgeI]
                 << ":: " << edges_[pEdges[edgeI]] << nl
                 << " edgeFaces:: " << edgeFaces_[pEdges[edgeI]] << nl
                 << " edgePoints:: " << edgePoints_[pEdges[edgeI]]
                 << endl;
        }

        Info << "Added faces: " << endl;
        forAll(newFaceIndex, faceI)
        {
            if (newFaceIndex[faceI] == -1)
            {
                continue;
            }

            Info << newFaceIndex[faceI] << ":: "
                 << faces_[newFaceIndex[faceI]]
                 << endl;
        }

        Info << "Added cells: " << endl;
        forAll(newCellIndex, cellI)
        {
            if (newCellIndex[cellI] == -1)
            {
                continue;
            }

            Info << newCellIndex[cellI] << ":: "
                 << cells_[newCellIndex[cellI]]
                 << endl;
        }

        // Write out VTK files after change
        if (debug > 3)
        {
            labelList vtkCells;

            if (cellsForRemoval[1] == -1)
            {
                vtkCells.setSize(3);

                // Fill in cell indices
                vtkCells[0] = newCellIndex[0];
                vtkCells[1] = newCellIndex[1];
                vtkCells[2] = newCellIndex[2];
            }
            else
            {
                vtkCells.setSize(6);

                // Fill in cell indices
                forAll(newCellIndex, indexI)
                {
                    vtkCells[indexI] = newCellIndex[indexI];
                }
            }

            writeVTK
            (
                Foam::name(fIndex)
              + "Trisect_1",
                vtkCells
            );
        }
    }

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    statistics_[3]++;

    // Increment the number of modifications
    statistics_[0]++;

    // Specify that the operation was successful
    map.type() = 1;

    // Return the changeMap
    return map;
}


// Utility method to compute the quality of a
// vertex hull around an edge after bisection.
scalar dynamicTopoFvMesh::computeBisectionQuality
(
    const label eIndex
) const
{
    scalar minQuality = GREAT, minVolume = GREAT;
    scalar cQuality = 0.0, oldVolume = 0.0;

    // Obtain a reference to this edge and corresponding edgePoints
    const edge& edgeToCheck = edges_[eIndex];
    const labelList& hullVertices = edgePoints_[eIndex];

    // Obtain point references
    const point& a = points_[edgeToCheck[0]];
    const point& c = points_[edgeToCheck[1]];

    const point& aOld = oldPoints_[edgeToCheck[0]];
    const point& cOld = oldPoints_[edgeToCheck[1]];

    // Compute the mid-point of the edge
    point midPoint = 0.5*(a + c);
    point oldPoint = 0.5*(aOld + cOld);

    if (whichEdgePatch(eIndex) < 0)
    {
        // Internal edge.
        forAll(hullVertices, indexI)
        {
            label prevIndex = hullVertices.rcIndex(indexI);

            // Pick vertices off the list
            const point& b = points_[hullVertices[prevIndex]];
            const point& d = points_[hullVertices[indexI]];

            const point& bOld = oldPoints_[hullVertices[prevIndex]];
            const point& dOld = oldPoints_[hullVertices[indexI]];

            // Compute the quality of the upper half.
            cQuality = tetMetric_(a, b, midPoint, d);

            // Compute old volume of the upper half.
            oldVolume = tetPointRef(aOld, bOld, oldPoint, dOld).mag();

            // Check if the volume / quality is worse
            minQuality = Foam::min(cQuality, minQuality);
            minVolume = Foam::min(oldVolume, minVolume);

            // Compute the quality of the lower half.
            cQuality = tetMetric_(midPoint, b, c, d);

            // Compute old volume of the lower half.
            oldVolume = tetPointRef(oldPoint, bOld, cOld, dOld).mag();

            // Check if the quality is worse
            minQuality = Foam::min(cQuality, minQuality);
            minVolume = Foam::min(oldVolume, minVolume);
        }
    }
    else
    {
        // Boundary edge.
        for(label indexI = 1; indexI < hullVertices.size(); indexI++)
        {
            // Pick vertices off the list
            const point& b = points_[hullVertices[indexI-1]];
            const point& d = points_[hullVertices[indexI]];

            const point& bOld = oldPoints_[hullVertices[indexI-1]];
            const point& dOld = oldPoints_[hullVertices[indexI]];

            // Compute the quality of the upper half.
            cQuality = tetMetric_(a, b, midPoint, d);

            // Compute old volume of the upper half.
            oldVolume = tetPointRef(aOld, bOld, oldPoint, dOld).mag();

            // Check if the quality is worse
            minQuality = Foam::min(cQuality, minQuality);
            minVolume = Foam::min(oldVolume, minVolume);

            // Compute the quality of the lower half.
            cQuality = tetMetric_(midPoint, b, c, d);

            // Compute old volume of the lower half.
            oldVolume = tetPointRef(oldPoint, bOld, cOld, dOld).mag();

            // Check if the quality is worse
            minQuality = Foam::min(cQuality, minQuality);
            minVolume = Foam::min(oldVolume, minVolume);
        }
    }

    // Ensure that the mesh is valid
    if (minQuality < sliverThreshold_)
    {
        if (debug > 3 && minQuality < 0.0)
        {
            // Write out cells for post processing.
            labelHashSet iCells;

            const labelList& eFaces = edgeFaces_[eIndex];

            forAll(eFaces, faceI)
            {
                if (!iCells.found(owner_[eFaces[faceI]]))
                {
                    iCells.insert(owner_[eFaces[faceI]]);
                }

                if (!iCells.found(neighbour_[eFaces[faceI]]))
                {
                    iCells.insert(neighbour_[eFaces[faceI]]);
                }
            }

            writeVTK(Foam::name(eIndex) + "_iCells", iCells.toc());
        }

        if (debug > 2)
        {
            InfoIn
            (
                "scalar dynamicTopoFvMesh::computeBisectionQuality"
                "(const label eIndex) const"
            )
                << "Bisecting edge will fall below the "
                << "sliver threshold of: " << sliverThreshold_ << nl
                << "Edge: " << eIndex << ": " << edgeToCheck << nl
                << "EdgePoints: " << hullVertices << nl
                << "Minimum Quality: " << minQuality << nl
                << "Mid point: " << midPoint
                << endl;
        }
    }

    // If a negative old-volume was encountered,
    // return an invalid quality.
    if (minVolume < 0.0)
    {
        return minVolume;
    }

    return minQuality;
}


// Utility method to compute the quality of cells
// around a face after trisection.
scalar dynamicTopoFvMesh::computeTrisectionQuality
(
    const label fIndex
) const
{
    scalar minQuality = GREAT;
    scalar cQuality = 0.0;

    point midPoint;

    // Fetch the midPoint
    midPoint = faces_[fIndex].centre(points_);

    FixedList<label,2> apexPoint(-1);

    // Find the apex point
    apexPoint[0] =
    (
        meshOps::tetApexPoint
        (
            owner_[fIndex],
            fIndex,
            faces_,
            cells_
        )
    );

    const face& faceToCheck = faces_[fIndex];

    forAll(faceToCheck, pointI)
    {
        // Pick vertices off the list
        const point& b = points_[faceToCheck[pointI]];
        const point& c = points_[apexPoint[0]];
        const point& d = points_[faceToCheck[faceToCheck.fcIndex(pointI)]];

        // Compute the quality of the upper half.
        cQuality = tetMetric_(midPoint, b, c, d);

        // Check if the quality is worse
        minQuality = Foam::min(cQuality, minQuality);
    }

    if (whichPatch(fIndex) == -1)
    {
        apexPoint[1] =
        (
            meshOps::tetApexPoint
            (
                neighbour_[fIndex],
                fIndex,
                faces_,
                cells_
            )
        );

        forAll(faceToCheck, pointI)
        {
            // Pick vertices off the list
            const point& b = points_[faceToCheck[pointI]];
            const point& c = points_[apexPoint[1]];
            const point& d = points_[faceToCheck[faceToCheck.rcIndex(pointI)]];

            // Compute the quality of the upper half.
            cQuality = tetMetric_(midPoint, b, c, d);

            // Check if the quality is worse
            minQuality = Foam::min(cQuality, minQuality);
        }
    }

    return minQuality;
}


// Split a set of internal faces into boundary faces
//   - Add boundary faces and edges to the patch specified by 'patchIndex'
//   - Cell color should specify a binary value dictating either side
//     of the split face.
void dynamicTopoFvMesh::splitInternalFaces
(
    const label patchIndex,
    const labelList& internalFaces,
    const Map<bool>& cellColors
)
{
    Map<label> mirrorPointLabels;
    FixedList<Map<label>, 2> mirrorEdgeLabels, mirrorFaceLabels;

    // First loop through the list and accumulate a list of
    // points and edges that need to be duplicated.
    forAll(internalFaces, faceI)
    {
        const face& faceToCheck = faces_[internalFaces[faceI]];

        forAll(faceToCheck, pointI)
        {
            if (!mirrorPointLabels.found(faceToCheck[pointI]))
            {
                mirrorPointLabels.insert(faceToCheck[pointI], -1);
            }
        }

        const labelList& fEdges = faceEdges_[internalFaces[faceI]];

        forAll(fEdges, edgeI)
        {
            if (!mirrorEdgeLabels[0].found(fEdges[edgeI]))
            {
                mirrorEdgeLabels[0].insert(fEdges[edgeI], -1);
            }
        }
    }

    // Now for every point in the list, add a new one.
    // Add a mapping entry as well.
    forAllIter(Map<label>, mirrorPointLabels, pIter)
    {
        // Obtain a copy of the point before adding it,
        // since the reference might become invalid during list resizing.
        point newPoint = points_[pIter.key()];
        point oldPoint = oldPoints_[pIter.key()];

        pIter() = insertPoint(newPoint, oldPoint, labelList(1, pIter.key()));

        if (!twoDMesh_)
        {
            const labelList& pEdges = pointEdges_[pIter.key()];

            labelHashSet edgesToRemove;

            forAll(pEdges, edgeI)
            {
                const labelList& eFaces = edgeFaces_[pEdges[edgeI]];

                bool allTrue = true;

                forAll(eFaces, faceI)
                {
                    label own = owner_[eFaces[faceI]];
                    label nei = neighbour_[eFaces[faceI]];

                    // Check if an owner/neighbour cell is false
                    if (!cellColors[own])
                    {
                        allTrue = false;
                        break;
                    }

                    if (nei != -1)
                    {
                        if (!cellColors[nei])
                        {
                            allTrue = false;
                            break;
                        }
                    }
                }

                if (allTrue)
                {
                    // Mark this edge label to be discarded later
                    edgesToRemove.insert(pEdges[edgeI]);
                }
            }

            // It is dangerous to use the pointEdges references,
            // so call it using array-lookup instead.
            forAllIter(labelHashSet, edgesToRemove, hsIter)
            {
                // Add the edge to the mirror point list
                meshOps::sizeUpList
                (
                    hsIter.key(),
                    pointEdges_[pIter()]
                );

                // Remove the edge from the original point list
                meshOps::sizeDownList
                (
                    hsIter.key(),
                    pointEdges_[pIter.key()]
                );
            }
        }
    }

    if (debug > 3)
    {
        label i = 0;
        labelList mPoints(mirrorPointLabels.size());

        if (!twoDMesh_)
        {
            forAllIter(Map<label>, mirrorPointLabels, pIter)
            {
                writeVTK
                (
                    "pEdges_o_" + Foam::name(pIter.key()) + '_',
                    pointEdges_[pIter.key()],
                    1
                );

                writeVTK
                (
                    "pEdges_m_" + Foam::name(pIter()) + '_',
                    pointEdges_[pIter()],
                    1
                );

                mPoints[i++] = pIter();
            }

            writeVTK
            (
                "points_o_",
                mirrorPointLabels.toc(),
                0
            );

            writeVTK
            (
                "points_m_",
                mPoints,
                0
            );
        }
    }

    // For every internal face, add a new one.
    //  - Stick to the rule:
    //    [1] Every cell marked false keeps the existing entities.
    //    [2] Every cell marked true gets new points/edges/faces.
    //  - If faces are improperly oriented, reverse them.
    forAll(internalFaces, faceI)
    {
        FixedList<face, 2> newFace;
        FixedList<label, 2> newFaceIndex(-1);
        FixedList<label, 2> newOwner(-1);

        label oldOwn = owner_[internalFaces[faceI]];
        label oldNei = neighbour_[internalFaces[faceI]];

        if (cellColors[oldOwn] && !cellColors[oldNei])
        {
            // The owner gets a new boundary face.
            // Note that orientation is already correct.
            newFace[0] = faces_[internalFaces[faceI]];

            // The neighbour needs to have its face reversed
            // and moved to the boundary patch, thereby getting
            // deleted in the process.
            newFace[1] = newFace[0].reverseFace();

            newOwner[0] = oldOwn;
            newOwner[1] = oldNei;
        }
        else
        if (!cellColors[oldOwn] && cellColors[oldNei])
        {
            // The neighbour gets a new boundary face.
            // The face is oriented in the opposite sense, however.
            newFace[0] = faces_[internalFaces[faceI]].reverseFace();

            // The owner keeps the existing face and orientation.
            // But it also needs to be moved to the boundary.
            newFace[1] = faces_[internalFaces[faceI]];

            newOwner[0] = oldNei;
            newOwner[1] = oldOwn;
        }
        else
        {
            // Something's wrong here.
            FatalErrorIn
            (
                "dynamicTopoFvMesh::splitInternalFaces\n"
                "(\n"
                "    const label patchIndex,\n"
                "    const labelList& internalFaces,\n"
                "    const Map<bool>& cellColors\n"
                ")\n"
            )
                << nl << " Face: "
                << internalFaces[faceI]
                << " has cells which are improperly marked: " << nl
                << oldOwn << ":: " << cellColors[oldOwn] << nl
                << oldNei << ":: " << cellColors[oldNei]
                << abort(FatalError);
        }

        // Renumber point labels for the first new face.
        forAll(newFace[0], pointI)
        {
            newFace[0][pointI] = mirrorPointLabels[newFace[0][pointI]];
        }

        // Insert the new boundary faces.
        forAll(newFace, indexI)
        {
            newFaceIndex[indexI] =
            (
                insertFace
                (
                    patchIndex,
                    newFace[indexI],
                    newOwner[indexI],
                    -1
                )
            );

            // Make an identical faceEdges entry.
            // This will be renumbered once new edges are added.
            labelList newFaceEdges(faceEdges_[internalFaces[faceI]]);

            faceEdges_.append(newFaceEdges);

            // Replace face labels on cells
            meshOps::replaceLabel
            (
                internalFaces[faceI],
                newFaceIndex[indexI],
                cells_[newOwner[indexI]]
            );

            // Make face mapping entries for posterity.
            mirrorFaceLabels[indexI].insert
            (
                internalFaces[faceI],
                newFaceIndex[indexI]
            );
        }
    }

    // For every edge in the list, add a new one.
    // We'll deal with correcting edgeFaces and edgePoints later.
    forAllIter(Map<label>, mirrorEdgeLabels[0], eIter)
    {
        // Obtain copies for the append method
        edge origEdge = edges_[eIter.key()];
        labelList newEdgeFaces(edgeFaces_[eIter.key()]);

        eIter() =
        (
            insertEdge
            (
                patchIndex,
                edge
                (
                    mirrorPointLabels[origEdge[0]],
                    mirrorPointLabels[origEdge[1]]
                ),
                newEdgeFaces,
                labelList(0)
            )
        );

        // Is the original edge an internal one?
        // If it is, we need to move it to the boundary.
        if (whichEdgePatch(eIter.key()) == -1)
        {
            label rplEdgeIndex =
            (
                insertEdge
                (
                    patchIndex,
                    origEdge,
                    newEdgeFaces,
                    labelList(0)
                )
            );

            // Map the new entry.
            mirrorEdgeLabels[1].insert(eIter.key(), rplEdgeIndex);
        }
        else
        {
            // This is already a boundary edge.
            // Make an identical map.
            mirrorEdgeLabels[1].insert(eIter.key(), eIter.key());
        }
    }

    // Correct edgeFaces for all new edges
    forAll(mirrorEdgeLabels, indexI)
    {
        forAllIter(Map<label>, mirrorEdgeLabels[indexI], eIter)
        {
            labelList& eFaces = edgeFaces_[eIter()];

            labelHashSet facesToRemove;

            forAll(eFaces, faceI)
            {
                bool remove = false;

                label own = owner_[eFaces[faceI]];
                label nei = neighbour_[eFaces[faceI]];

                if
                (
                    (!cellColors[own] && !indexI) ||
                    ( cellColors[own] &&  indexI)
                )
                {
                    remove = true;
                }

                if (nei != -1)
                {
                    if
                    (
                        (!cellColors[nei] && !indexI) ||
                        ( cellColors[nei] &&  indexI)
                    )
                    {
                        remove = true;
                    }
                }

                if (mirrorFaceLabels[indexI].found(eFaces[faceI]))
                {
                    // Perform a replacement instead of a removal.
                    eFaces[faceI] = mirrorFaceLabels[indexI][eFaces[faceI]];

                    remove = false;
                }

                if (remove)
                {
                    facesToRemove.insert(eFaces[faceI]);
                }
            }

            // Now successively size down edgeFaces.
            // It is dangerous to use the eFaces reference,
            // so call it using array-lookup.
            forAllIter(labelHashSet, facesToRemove, hsIter)
            {
                meshOps::sizeDownList(hsIter.key(), edgeFaces_[eIter()]);
            }
        }
    }

    // Renumber faceEdges for all faces connected to new edges
    forAll(mirrorEdgeLabels, indexI)
    {
        forAllIter(Map<label>, mirrorEdgeLabels[indexI], eIter)
        {
            const labelList& eFaces = edgeFaces_[eIter()];

            forAll(eFaces, faceI)
            {
                labelList& fEdges = faceEdges_[eFaces[faceI]];

                forAll(fEdges, edgeI)
                {
                    if (mirrorEdgeLabels[indexI].found(fEdges[edgeI]))
                    {
                        fEdges[edgeI] =
                        (
                            mirrorEdgeLabels[indexI][fEdges[edgeI]]
                        );
                    }
                }
            }
        }
    }

    if (twoDMesh_)
    {
        // Renumber edges and faces
        forAllIter(Map<label>, mirrorEdgeLabels[0], eIter)
        {
            const labelList& eFaces = edgeFaces_[eIter()];

            // Two levels of indirection to ensure
            // that all entities we renumbered.
            // A flip-side for the lack of a pointEdges list in 2D.
            forAll(eFaces, faceI)
            {
                const labelList& fEdges = faceEdges_[eFaces[faceI]];

                forAll(fEdges, edgeI)
                {
                    // Renumber this edge.
                    edge& edgeToCheck = edges_[fEdges[edgeI]];

                    forAll(edgeToCheck, pointI)
                    {
                        if (mirrorPointLabels.found(edgeToCheck[pointI]))
                        {
                            edgeToCheck[pointI] =
                            (
                                mirrorPointLabels[edgeToCheck[pointI]]
                            );
                        }
                    }

                    // Also renumber faces connected to this edge.
                    const labelList& efFaces = edgeFaces_[fEdges[edgeI]];

                    forAll(efFaces, faceJ)
                    {
                        face& faceToCheck = faces_[efFaces[faceJ]];

                        forAll(faceToCheck, pointI)
                        {
                            if
                            (
                                mirrorPointLabels.found(faceToCheck[pointI])
                            )
                            {
                                faceToCheck[pointI] =
                                (
                                    mirrorPointLabels[faceToCheck[pointI]]
                                );
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        // Point renumbering of entities connected to mirror points
        forAllIter(Map<label>, mirrorPointLabels, pIter)
        {
            const labelList& pEdges = pointEdges_[pIter()];

            forAll(pEdges, edgeI)
            {
                // Renumber this edge.
                edge& edgeToCheck = edges_[pEdges[edgeI]];

                forAll(edgeToCheck, pointI)
                {
                    if (edgeToCheck[pointI] == pIter.key())
                    {
                        edgeToCheck[pointI] = pIter();
                    }
                }

                // Also renumber faces connected to this edge.
                const labelList& eFaces = edgeFaces_[pEdges[edgeI]];

                forAll(eFaces, faceI)
                {
                    face& faceToCheck = faces_[eFaces[faceI]];

                    forAll(faceToCheck, pointI)
                    {
                        if (faceToCheck[pointI] == pIter.key())
                        {
                            faceToCheck[pointI] = pIter();
                        }
                    }
                }
            }
        }

        // Scan edges connected to mirror points,
        // and correct any edgePoints entries for
        // edges connected to the other vertex.
        forAllIter(Map<label>, mirrorPointLabels, pIter)
        {
            const labelList& pEdges = pointEdges_[pIter()];

            forAll(pEdges, edgeI)
            {
                label otherVertex = edges_[pEdges[edgeI]].otherVertex(pIter());

                // Scan edgePoints for edges connected to this point
                const labelList& opEdges = pointEdges_[otherVertex];

                forAll(opEdges, edgeJ)
                {
                    labelList& ePoints = edgePoints_[opEdges[edgeJ]];

                    forAll(ePoints, pointI)
                    {
                        if (mirrorPointLabels.found(ePoints[pointI]))
                        {
                            // Replace this point with the mirror point
                            ePoints[pointI] =
                            (
                                mirrorPointLabels[ePoints[pointI]]
                            );
                        }
                    }
                }
            }
        }

        // Build edgePoints for new edges
        forAll(mirrorEdgeLabels, indexI)
        {
            forAllIter(Map<label>, mirrorEdgeLabels[indexI], eIter)
            {
                buildEdgePoints(eIter());
            }
        }
    }

    // Now that we're done with the internal faces, remove them.
    forAll(internalFaces, faceI)
    {
        removeFace(internalFaces[faceI]);
    }

    // Remove old internal edges as well.
    forAllIter(Map<label>, mirrorEdgeLabels[1], eIter)
    {
        if (eIter.key() != eIter())
        {
            removeEdge(eIter.key());
        }
    }

    // Set the flag
    topoChangeFlag_ = true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
