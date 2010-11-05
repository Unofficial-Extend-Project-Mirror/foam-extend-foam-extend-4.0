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
#include "triPointRef.H"
#include "linePointRef.H"
#include "multiThreader.H"
#include "dynamicTopoFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Perform a Delaunay test on an internal face
bool dynamicTopoFvMesh::testDelaunay
(
    const label fIndex
) const
{
    bool failed = false;
    label eIndex = -1, pIndex = -1, fLabel = -1;
    FixedList<bool,2> foundTriFace(false);
    FixedList<FixedList<label,3>,2> triFaces(FixedList<label,3>(-1));

    // Boundary faces are discarded.
    if (whichPatch(fIndex) > -1)
    {
        {
            return failed;
        }
    }

    {
        const labelList& fEdges = faceEdges_[fIndex];

        forAll(fEdges, edgeI)
        {
            // Break out if both triangular faces are found
            if (foundTriFace[0] && foundTriFace[1])
            {
                break;
            }

            // Obtain edgeFaces for this edge
            const labelList& eFaces = edgeFaces_[fEdges[edgeI]];

            forAll(eFaces, faceI)
            {
                const face& thisFace = faces_[eFaces[faceI]];

                if (thisFace.size() == 3)
                {
                    if (foundTriFace[0])
                    {
                        // Update the second face.
                        triFaces[1][0] = thisFace[0];
                        triFaces[1][1] = thisFace[1];
                        triFaces[1][2] = thisFace[2];

                        foundTriFace[1] = true;

                        // Take this edge
                        eIndex = fEdges[edgeI];
                    }
                    else
                    {
                        // Update the first face.
                        triFaces[0][0] = thisFace[0];
                        triFaces[0][1] = thisFace[1];
                        triFaces[0][2] = thisFace[2];

                        foundTriFace[0] = true;

                        fLabel = eFaces[faceI];
                    }
                }
            }
        }
    }

    // Obtain point references for the first face
    point a = points_[triFaces[0][0]];

    const face& faceToCheck = faces_[fLabel];

    vector cCentre =
    (
        triPointRef
        (
            points_[faceToCheck[0]],
            points_[faceToCheck[1]],
            points_[faceToCheck[2]]
        ).circumCentre()
    );

    point otherPoint = vector::zero;
    scalar rSquared = (a - cCentre)&(a - cCentre);

    // Find the isolated point on the second face
    {
        const edge& e = edges_[eIndex];

        // Check the first point
        if (triFaces[1][0] != e.start() && triFaces[1][0] != e.end())
        {
            pIndex = triFaces[1][0];
        }

        // Check the second point
        if (triFaces[1][1] != e.start() && triFaces[1][1] != e.end())
        {
            pIndex = triFaces[1][1];
        }

        // Check the third point
        if (triFaces[1][2] != e.start() && triFaces[1][2] != e.end())
        {
            pIndex = triFaces[1][2];
        }

        // ...and determine whether it lies in this circle
        otherPoint = points_[pIndex];
    }

    if (((otherPoint - cCentre)&(otherPoint - cCentre)) < rSquared)
    {
        // Failed the test.
        failed = true;
    }

    return failed;
}


// Method for the swapping of a quad-face in 2D
// - Returns a changeMap with a type specifying:
//     1: Swap sequence was successful
//    -1: Swap sequence failed
const changeMap dynamicTopoFvMesh::swapQuadFace
(
    const label fIndex
)
{
    changeMap map;

    face f;
    bool found = false;
    label commonIndex = -1;
    FixedList<label,4> otherPointIndex(-1), nextToOtherPoint(-1);
    FixedList<label,2> c0BdyIndex, c0IntIndex, c1BdyIndex, c1IntIndex;
    FixedList<face,2>  c0BdyFace,  c0IntFace,  c1BdyFace,  c1IntFace;
    FixedList<label,2> commonEdgeIndex(-1);
    FixedList<edge,2>  commonEdges;
    FixedList<label,4> otherEdgeIndex(-1);
    FixedList<label,4> commonFaceIndex(-1), cornerEdgeIndex(-1);
    FixedList<face,4>  commonFaces(face(3)), commonIntFaces(face(4));
    FixedList<label,4> commonIntFaceIndex(-1);
    FixedList<bool,2> foundTriFace0(false), foundTriFace1(false);
    FixedList<face,2> triFaces0(face(3)), triFaces1(face(3));

    // Get the two cells on either side...
    label c0 = owner_[fIndex];
    label c1 = neighbour_[fIndex];

    // Prepare two new cells
    FixedList<cell, 2> newCell(cell(5));

    // Need to find common faces and edges...
    // At the end of this loop, commonFaces [0] & [1] share commonEdge[0]
    // and commonFaces [2] & [3] share commonEdge[1]
    // Also, commonFaces[0] & [2] are connected to cell[0],
    // and commonFaces[1] & [3] are connected to cell[1]
    const labelList& fEdges = faceEdges_[fIndex];

    forAll(fEdges, edgeI)
    {
        // Break out if all triangular faces are found
        if
        (
            foundTriFace0[0] && foundTriFace0[1] &&
            foundTriFace1[0] && foundTriFace1[1]
        )
        {
            break;
        }

        // Obtain edgeFaces for this edge
        const labelList& eFaces = edgeFaces_[fEdges[edgeI]];

        forAll(eFaces, faceI)
        {
            const face& eFace = faces_[eFaces[faceI]];

            if (eFace.size() == 3)
            {
                // Found a triangular face. Determine which cell it belongs to.
                if (owner_[eFaces[faceI]] == c0)
                {
                    if (foundTriFace0[0])
                    {
                        // Update the second face on cell[0].
                        commonIndex = 2;
                        foundTriFace0[1] = true;

                        if (foundTriFace1[1])
                        {
                            commonEdgeIndex[1] = fEdges[edgeI];
                            commonEdges[1] = edges_[fEdges[edgeI]];
                        }
                    }
                    else
                    {
                        // Update the first face on cell[0].
                        commonIndex = 0;
                        foundTriFace0[0] = true;

                        if (foundTriFace1[0])
                        {
                            commonEdgeIndex[0] = fEdges[edgeI];
                            commonEdges[0] = edges_[fEdges[edgeI]];
                        }
                    }
                }
                else
                {
                    if (foundTriFace1[0])
                    {
                        // Update the second face on cell[1].
                        commonIndex = 3;
                        foundTriFace1[1] = true;

                        if (foundTriFace0[1])
                        {
                            commonEdgeIndex[1] = fEdges[edgeI];
                            commonEdges[1] = edges_[fEdges[edgeI]];
                        }
                    }
                    else
                    {
                        // Update the first face on cell[1].
                        commonIndex = 1;
                        foundTriFace1[0] = true;

                        if (foundTriFace0[0])
                        {
                            commonEdgeIndex[0] = fEdges[edgeI];
                            commonEdges[0] = edges_[fEdges[edgeI]];
                        }
                    }
                }

                // Store the face and index
                commonFaces[commonIndex][0] = eFace[0];
                commonFaces[commonIndex][1] = eFace[1];
                commonFaces[commonIndex][2] = eFace[2];

                commonFaceIndex[commonIndex] = eFaces[faceI];
            }
        }
    }

    // Find the interior/boundary faces.
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

    // Find the points that don't lie on shared edges
    // and the points next to them (for orientation)
    meshOps::findIsolatedPoint
    (
        commonFaces[1],
        commonEdges[0],
        otherPointIndex[1],
        nextToOtherPoint[1]
    );

    meshOps::findIsolatedPoint
    (
        commonFaces[0],
        commonEdges[0],
        otherPointIndex[0],
        nextToOtherPoint[0]
    );

    meshOps::findIsolatedPoint
    (
        commonFaces[2],
        commonEdges[1],
        otherPointIndex[2],
        nextToOtherPoint[2]
    );

    meshOps::findIsolatedPoint
    (
        commonFaces[3],
        commonEdges[1],
        otherPointIndex[3],
        nextToOtherPoint[3]
    );

    // Ensure that the configuration will be valid
    // using old point-positions. A simple area-based
    // calculation should suffice.
    FixedList<face, 2> triFaceOldPoints(face(3));

    triFaceOldPoints[0][0] = otherPointIndex[0];
    triFaceOldPoints[0][1] = nextToOtherPoint[0];
    triFaceOldPoints[0][2] = otherPointIndex[1];

    triFaceOldPoints[1][0] = otherPointIndex[1];
    triFaceOldPoints[1][1] = nextToOtherPoint[1];
    triFaceOldPoints[1][2] = otherPointIndex[0];

    // Assume XY plane here
    vector n = vector(0,0,1);

    forAll(triFaceOldPoints, faceI)
    {
        // Assume that centre-plane passes through the origin.
        vector xf, nf;

        xf = triFaceOldPoints[faceI].centre(oldPoints_);
        nf = triFaceOldPoints[faceI].normal(oldPoints_);

        if ((((xf & n) * n) & nf) < 0.0)
        {
            // This will yield an inverted cell. Bail out.
            return map;
        }
    }

    // Find the other two edges on the face being flipped
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

    // At the end of this loop, commonIntFaces [0] & [1] share otherEdges[0]
    // and commonIntFaces [2] & [3] share the otherEdges[1],
    // where [0],[2] lie on cell[0] and [1],[3] lie on cell[1]
    found = false;

    const labelList& e1 = faceEdges_[c0IntIndex[0]];

    forAll(e1,edgeI)
    {
        if (e1[edgeI] == otherEdgeIndex[0])
        {
            commonIntFaces[0] = c0IntFace[0];
            commonIntFaces[2] = c0IntFace[1];
            commonIntFaceIndex[0] = c0IntIndex[0];
            commonIntFaceIndex[2] = c0IntIndex[1];
            found = true; break;
        }
    }

    if (!found)
    {
        // The edge was not found before
        commonIntFaces[0] = c0IntFace[1];
        commonIntFaces[2] = c0IntFace[0];
        commonIntFaceIndex[0] = c0IntIndex[1];
        commonIntFaceIndex[2] = c0IntIndex[0];
    }

    found = false;

    const labelList& e3 = faceEdges_[c1IntIndex[0]];

    forAll(e3,edgeI)
    {
        if (e3[edgeI] == otherEdgeIndex[0])
        {
            commonIntFaces[1] = c1IntFace[0];
            commonIntFaces[3] = c1IntFace[1];
            commonIntFaceIndex[1] = c1IntIndex[0];
            commonIntFaceIndex[3] = c1IntIndex[1];
            found = true; break;
        }
    }

    if (!found)
    {
        // The edge was not found before
        commonIntFaces[1] = c1IntFace[1];
        commonIntFaces[3] = c1IntFace[0];
        commonIntFaceIndex[1] = c1IntIndex[1];
        commonIntFaceIndex[3] = c1IntIndex[0];
    }

    // Find two common edges between quad/quad faces
    meshOps::findCommonEdge
    (
        c0IntIndex[0],
        c0IntIndex[1],
        faceEdges_,
        otherEdgeIndex[2]
    );

    meshOps::findCommonEdge
    (
        c1IntIndex[0],
        c1IntIndex[1],
        faceEdges_,
        otherEdgeIndex[3]
    );

    // Find four common edges between quad/tri faces
    meshOps::findCommonEdge
    (
        commonFaceIndex[1],
        commonIntFaceIndex[1],
        faceEdges_,
        cornerEdgeIndex[0]
    );

    meshOps::findCommonEdge
    (
        commonFaceIndex[3],
        commonIntFaceIndex[1],
        faceEdges_,
        cornerEdgeIndex[1]
    );

    meshOps::findCommonEdge
    (
        commonFaceIndex[0],
        commonIntFaceIndex[2],
        faceEdges_,
        cornerEdgeIndex[2]
    );

    meshOps::findCommonEdge
    (
        commonFaceIndex[2],
        commonIntFaceIndex[2],
        faceEdges_,
        cornerEdgeIndex[3]
    );

    // Perform a few debug calls before modifications
    if (debug > 1)
    {
        Info << nl << nl << "Face: " << fIndex
             << " needs to be flipped. " << endl;

        Info << "Cell[0]: " << c0 << ": " << cells_[c0] << endl;
        Info << "Cell[1]: " << c1 << ": " << cells_[c1] << endl;

        if (debug > 2)
        {
            Info << "Common Faces: Set 1: "
                 << commonFaceIndex[0] << ": " << commonFaces[0] << ", "
                 << commonFaceIndex[1] << ": " << commonFaces[1] << endl;

            Info << "Common Faces: Set 2: "
                 << commonFaceIndex[2] << ": " << commonFaces[2] << ", "
                 << commonFaceIndex[3] << ": " << commonFaces[3] << endl;

            Info << "Old face: " << faces_[fIndex] << endl;
            Info << "Old faceEdges: " << faceEdges_[fIndex] << endl;
        }

        // Write out VTK files before change
        if (debug > 3)
        {
            labelList cellHull(2, -1);

            cellHull[0] = c0;
            cellHull[1] = c1;

            writeVTK
            (
                Foam::name(fIndex)
              + "_Swap_0",
                cellHull
            );
        }
    }

    // Modify the five faces belonging to this hull
    face newFace = faces_[fIndex];
    labelList newFEdges = faceEdges_[fIndex];
    FixedList<face, 4> newBdyFace(face(3));
    FixedList<edge, 2> newEdges;

    // Assign current faces
    forAll(newBdyFace, faceI)
    {
        newBdyFace[faceI] = faces_[commonFaceIndex[faceI]];
    }

    label c0count=0, c1count=0;

    // Size down edgeFaces for the original face
    meshOps::sizeDownList
    (
        fIndex,
        edgeFaces_[otherEdgeIndex[0]]
    );

    meshOps::sizeDownList
    (
        fIndex,
        edgeFaces_[otherEdgeIndex[1]]
    );

    // Size up edgeFaces for the face after flipping
    meshOps::sizeUpList
    (
        fIndex,
        edgeFaces_[otherEdgeIndex[2]]
    );

    meshOps::sizeUpList
    (
        fIndex,
        edgeFaces_[otherEdgeIndex[3]]
    );

    // Replace edgeFaces and faceEdges for the
    // four (out of 8 total) corner edges of this hull.
    meshOps::replaceLabel
    (
        cornerEdgeIndex[0],
        cornerEdgeIndex[2],
        faceEdges_[commonFaceIndex[1]]
    );

    meshOps::replaceLabel
    (
        commonFaceIndex[1],
        commonFaceIndex[0],
        edgeFaces_[cornerEdgeIndex[0]]
    );

    meshOps::replaceLabel
    (
        cornerEdgeIndex[1],
        cornerEdgeIndex[3],
        faceEdges_[commonFaceIndex[3]]
    );

    meshOps::replaceLabel
    (
        commonFaceIndex[3],
        commonFaceIndex[2],
        edgeFaces_[cornerEdgeIndex[1]]
    );

    meshOps::replaceLabel
    (
        cornerEdgeIndex[2],
        cornerEdgeIndex[0],
        faceEdges_[commonFaceIndex[0]]
    );

    meshOps::replaceLabel
    (
        commonFaceIndex[0],
        commonFaceIndex[1],
        edgeFaces_[cornerEdgeIndex[2]]
    );

    meshOps::replaceLabel
    (
        cornerEdgeIndex[3],
        cornerEdgeIndex[1],
        faceEdges_[commonFaceIndex[2]]
    );

    meshOps::replaceLabel
    (
        commonFaceIndex[2],
        commonFaceIndex[3],
        edgeFaces_[cornerEdgeIndex[3]]
    );

    // Define parameters for the new flipped face
    newFace[0] = otherPointIndex[0];
    newFace[1] = otherPointIndex[1];
    newFace[2] = otherPointIndex[3];
    newFace[3] = otherPointIndex[2];

    newFEdges[0] = otherEdgeIndex[2];
    newFEdges[1] = commonEdgeIndex[0];
    newFEdges[2] = otherEdgeIndex[3];
    newFEdges[3] = commonEdgeIndex[1];

    newCell[0][c0count++] = fIndex;
    newCell[1][c1count++] = fIndex;

    owner_[fIndex] = c0;
    neighbour_[fIndex] = c1;

    // Both commonEdges need to be renumbered.
    newEdges[0][0] = otherPointIndex[0];
    newEdges[0][1] = otherPointIndex[1];

    newEdges[1][0] = otherPointIndex[3];
    newEdges[1][1] = otherPointIndex[2];

    // Four modified boundary faces need to be constructed,
    // but right-handedness is also important.
    // Take a cue from the existing boundary-face orientation

    // Zeroth boundary face - Owner c[0], Neighbour -1
    newBdyFace[0][0] = otherPointIndex[0];
    newBdyFace[0][1] = nextToOtherPoint[0];
    newBdyFace[0][2] = otherPointIndex[1];
    newCell[0][c0count++] = commonFaceIndex[0];
    owner_[commonFaceIndex[0]] = c0;
    neighbour_[commonFaceIndex[0]] = -1;

    // First boundary face - Owner c[1], Neighbour -1
    newBdyFace[1][0] = otherPointIndex[1];
    newBdyFace[1][1] = nextToOtherPoint[1];
    newBdyFace[1][2] = otherPointIndex[0];
    newCell[1][c1count++] = commonFaceIndex[1];
    owner_[commonFaceIndex[1]] = c1;
    neighbour_[commonFaceIndex[1]] = -1;

    // Second boundary face - Owner c[0], Neighbour -1
    newBdyFace[2][0] = otherPointIndex[3];
    newBdyFace[2][1] = nextToOtherPoint[3];
    newBdyFace[2][2] = otherPointIndex[2];
    newCell[0][c0count++] = commonFaceIndex[2];
    owner_[commonFaceIndex[2]] = c0;
    neighbour_[commonFaceIndex[2]] = -1;

    // Third boundary face - Owner c[1], Neighbour -1
    newBdyFace[3][0] = otherPointIndex[2];
    newBdyFace[3][1] = nextToOtherPoint[2];
    newBdyFace[3][2] = otherPointIndex[3];
    newCell[1][c1count++] = commonFaceIndex[3];
    owner_[commonFaceIndex[3]] = c1;
    neighbour_[commonFaceIndex[3]] = -1;

    if (debug > 1)
    {
        Info << "New flipped face: " << newFace << endl;

        if (debug > 2)
        {
            forAll(newBdyFace, faceI)
            {
                Info << "New boundary face[" << faceI << "]: "
                     << commonFaceIndex[faceI]
                     << ": " << newBdyFace[faceI] << endl;
            }
        }
    }

    // Check the orientation of the two quad faces, and modify as necessary
    label newOwn=0, newNei=0;

    // The quad face belonging to cell[1] now becomes a part of cell[0]
    if (neighbour_[commonIntFaceIndex[1]] == -1)
    {
        // Boundary face
        // Face doesn't need to be flipped, just update the owner
        f          = commonIntFaces[1];
        newOwn     = c0;
        newNei     = -1;
    }
    else
    if (owner_[commonIntFaceIndex[1]] == c1)
    {
        // This face is on the interior, check for previous owner
        // Upper-triangular ordering has to be maintained, however...
        if (c0 > neighbour_[commonIntFaceIndex[1]])
        {
            // Flip is necessary
            f          = commonIntFaces[1].reverseFace();
            newOwn     = neighbour_[commonIntFaceIndex[1]];
            newNei     = c0;

            setFlip(commonIntFaceIndex[1]);
        }
        else
        {
            // Flip isn't necessary, just change the owner
            f          = commonIntFaces[1];
            newOwn     = c0;
            newNei     = neighbour_[commonIntFaceIndex[1]];
        }
    }
    else
    if (neighbour_[commonIntFaceIndex[1]] == c1)
    {
        // This face is on the interior, check for previous neighbour
        // Upper-triangular ordering has to be maintained, however...
        if (c0 < owner_[commonIntFaceIndex[1]])
        {
            // Flip is necessary
            f          = commonIntFaces[1].reverseFace();
            newOwn     = c0;
            newNei     = owner_[commonIntFaceIndex[1]];

            setFlip(commonIntFaceIndex[1]);
        }
        else
        {
            // Flip isn't necessary, just change the neighbour
            f          = commonIntFaces[1];
            newOwn     = owner_[commonIntFaceIndex[1]];
            newNei     = c0;
        }
    }

    faces_[commonIntFaceIndex[1]] = f;
    newCell[0][c0count++] = commonIntFaceIndex[0];
    newCell[0][c0count++] = commonIntFaceIndex[1];
    owner_[commonIntFaceIndex[1]] = newOwn;
    neighbour_[commonIntFaceIndex[1]] = newNei;

    // The quad face belonging to cell[0] now becomes a part of cell[1]
    if (neighbour_[commonIntFaceIndex[2]] == -1)
    {
        // Boundary face
        // Face doesn't need to be flipped, just update the owner
        f          = commonIntFaces[2];
        newOwn     = c1;
        newNei     = -1;
    }
    else
    if (owner_[commonIntFaceIndex[2]] == c0)
    {
        // This face is on the interior, check for previous owner
        // Upper-triangular ordering has to be maintained, however...
        if (c1 > neighbour_[commonIntFaceIndex[2]])
        {
            // Flip is necessary
            f          = commonIntFaces[2].reverseFace();
            newOwn     = neighbour_[commonIntFaceIndex[2]];
            newNei     = c1;

            setFlip(commonIntFaceIndex[2]);
        }
        else
        {
            // Flip isn't necessary, just change the owner
            f          = commonIntFaces[2];
            newOwn     = c1;
            newNei     = neighbour_[commonIntFaceIndex[2]];
        }
    }
    else
    if (neighbour_[commonIntFaceIndex[2]] == c0)
    {
        // This face is on the interior, check for previous neighbour
        // Upper-triangular ordering has to be maintained, however...
        if (c1 < owner_[commonIntFaceIndex[2]])
        {
            // Flip is necessary
            f          = commonIntFaces[2].reverseFace();
            newOwn     = c1;
            newNei     = owner_[commonIntFaceIndex[2]];

            setFlip(commonIntFaceIndex[2]);
        }
        else
        {
            // Flip isn't necessary, just change the neighbour
            f          = commonIntFaces[2];
            newOwn     = owner_[commonIntFaceIndex[2]];
            newNei     = c1;
        }
    }

    faces_[commonIntFaceIndex[2]] = f;
    newCell[1][c1count++] = commonIntFaceIndex[2];
    newCell[1][c1count++] = commonIntFaceIndex[3];
    owner_[commonIntFaceIndex[2]] = newOwn;
    neighbour_[commonIntFaceIndex[2]] = newNei;

    // Now that all entities are properly configured,
    // overwrite the entries in connectivity lists.
    cells_[c0] = newCell[0];
    cells_[c1] = newCell[1];

    faces_[fIndex] = newFace;
    faceEdges_[fIndex] = newFEdges;

    forAll(newBdyFace, faceI)
    {
        faces_[commonFaceIndex[faceI]] = newBdyFace[faceI];
    }

    edges_[commonEdgeIndex[0]] = newEdges[0];
    edges_[commonEdgeIndex[1]] = newEdges[1];

    // Fill-in mapping information
    labelList mC(2, -1);
    mC[0] = c0; mC[1] = c1;

    forAll(mC, cellI)
    {
        // Set the mapping for this cell
        setCellMapping(mC[cellI], mC);
    }

    // Interpolate new fluxes for the flipped face.
    setFaceMapping(fIndex);

    // Write out VTK files after change
    if (debug > 3)
    {
        labelList cellHull(2, -1);

        cellHull[0] = c0;
        cellHull[1] = c1;

        writeVTK
        (
            Foam::name(fIndex)
          + "_Swap_1",
            cellHull
        );
    }

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    statistics_[1]++;

    // Return a successful operation.
    map.type() = 1;

    return map;
}


// Method for the swapping of an edge in 3D
//  - To be used mainly for testing purposes only.
//  - Use swap3DEdges on the entire mesh for efficiency.
void dynamicTopoFvMesh::swapEdge
(
    const label eIndex,
    bool forceOp
)
{
    // Dynamic programming variables
    labelList m;
    PtrList<scalarListList> Q;
    PtrList<labelListList> K, triangulations;

    // Allocate dynamic programming tables
    initTables(m, Q, K, triangulations);

    // Compute the minimum quality of cells around this edge
    scalar minQuality = computeMinQuality(eIndex);

    // Check if this edge is on a bounding curve
    if (checkBoundingCurve(eIndex))
    {
        FatalErrorIn
        (
            "void dynamicTopoFvMesh::swapEdge"
            "(const label eIndex, bool forceOp)"
        )
            << nl << " Cannot swap edges on bounding curves. "
            << abort(FatalError);
    }

    // Fill the dynamic programming tables
    if (fillTables(eIndex, minQuality, m, Q, K, triangulations))
    {
        // Check if edge-swapping is required.
        scalar newQuality = Q[0][0][m[0]-1];

        if (newQuality > minQuality)
        {
            // Remove this edge according to the swap sequence
            removeEdgeFlips(eIndex, minQuality, K, triangulations);
        }
        else
        if (forceOp)
        {
            if (newQuality < 0.0)
            {
                FatalErrorIn
                (
                    "void dynamicTopoFvMesh::swapEdge"
                    "(const label eIndex, bool forceOp)"
                )
                    << " Forcing swap on edge: " << eIndex
                    << ":: " << edges_[eIndex]
                    << " will yield an invalid cell quality: "
                    << newQuality << " Old Quality: " << minQuality
                    << abort(FatalError);
            }
            else
            {
                removeEdgeFlips(eIndex, minQuality, K, triangulations);
            }
        }
    }
}


// Allocate dynamic programming tables
void dynamicTopoFvMesh::initTables
(
    labelList& m,
    PtrList<scalarListList>& Q,
    PtrList<labelListList>& K,
    PtrList<labelListList>& triangulations,
    const label checkIndex
) const
{
    label mMax = maxTetsPerEdge_;

    // Check if resizing is necessary only for a particular index.
    if (checkIndex != -1)
    {
        m[checkIndex] = -1;
        Q[checkIndex].setSize((mMax-2),scalarList(mMax,-1.0));
        K[checkIndex].setSize((mMax-2),labelList(mMax,-1));
        triangulations[checkIndex].setSize(3,labelList((mMax-2),-1));

        return;
    }

    // Size all elements by default.
    label numIndices = -1;

    {
        numIndices = 1;
    }

    m.setSize(numIndices, -1);
    Q.setSize(numIndices);
    K.setSize(numIndices);
    triangulations.setSize(numIndices);

    forAll(Q, indexI)
    {
        Q.set
        (
            indexI,
            new scalarListList((mMax-2),scalarList(mMax,-1.0))
        );

        K.set
        (
            indexI,
            new labelListList((mMax-2),labelList(mMax,-1))
        );

        triangulations.set
        (
            indexI,
            new labelListList(3,labelList((mMax-2),-1))
        );
    }
}


// Utility method to fill the dynamic programming tables
//  - Returns true if the operation completed successfully.
//  - Returns false if tables could not be resized.
bool dynamicTopoFvMesh::fillTables
(
    const label eIndex,
    const scalar minQuality,
    labelList& m,
    PtrList<scalarListList>& Q,
    PtrList<labelListList>& K,
    PtrList<labelListList>& triangulations,
    const label checkIndex
) const
{
    const edge& edgeToCheck = edges_[eIndex];
    const labelList& hullVertices = edgePoints_[eIndex];

    // Fill in the size
    m[checkIndex] = hullVertices.size();

    // Check if a table-resize is necessary
    if (m[checkIndex] > maxTetsPerEdge_)
    {
        if (allowTableResize_)
        {
            // Resize the tables to account for
            // more tets per edge
            label& mtpe = const_cast<label&>(maxTetsPerEdge_);

            mtpe = m[checkIndex];

            // Clear tables for this index.
            Q[checkIndex].clear();
            K[checkIndex].clear();
            triangulations[checkIndex].clear();

            // Resize for this index.
            initTables(m, Q, K, triangulations, checkIndex);
        }
        else
        {
            // Can't resize. Bail out.
            return false;
        }
    }

    for (label i = (m[checkIndex]-3); i >= 0; i--)
    {
        for (label j = i+2; j < m[checkIndex]; j++)
        {
            for (label k = i+1; k < j; k++)
            {
                scalar q = (*tetMetric_)
                (
                    points_[hullVertices[i]],
                    points_[hullVertices[k]],
                    points_[hullVertices[j]],
                    points_[edgeToCheck[0]]
                );

                // For efficiency, check the bottom triangulation
                // only when the top one if less than the hull quality.
                if (q > minQuality)
                {
                    q =
                    (
                        Foam::min
                        (
                            q,
                            (*tetMetric_)
                            (
                                points_[hullVertices[j]],
                                points_[hullVertices[k]],
                                points_[hullVertices[i]],
                                points_[edgeToCheck[1]]
                            )
                        )
                    );
                }

                if (k < j-1)
                {
                    q = Foam::min(q,Q[checkIndex][k][j]);
                }

                if (k > i+1)
                {
                    q = Foam::min(q,Q[checkIndex][i][k]);
                }

                if ((k == i+1) || (q > Q[checkIndex][i][j]))
                {
                    Q[checkIndex][i][j] = q;
                    K[checkIndex][i][j] = k;
                }
            }
        }
    }

    // Print out tables for debugging
    if (debug > 3)
    {
        printTables(m, Q, K, checkIndex);
    }

    return true;
}


// Remove the edge according to the swap sequence.
// - Returns a changeMap with a type specifying:
//     1: Swap sequence was successful
//    -1: Swap sequence failed
const changeMap dynamicTopoFvMesh::removeEdgeFlips
(
    const label eIndex,
    const scalar minQuality,
    const PtrList<labelListList>& K,
    PtrList<labelListList>& triangulations,
    const label checkIndex
)
{
    changeMap map, slaveMap;

    if (debug > 2)
    {
        Info << " Removing edge : " << eIndex << " by flipping."
             << " Edge: " << edges_[eIndex]
             << " minQuality: " << minQuality << endl;
    }

    // Make a copy of edgePoints, since it will be
    // modified during swaps
    labelList hullVertices(edgePoints_[eIndex]);

    label m = hullVertices.size();

    labelList hullFaces(m, -1);
    labelList hullCells(m, -1);
    labelList hullEdges(m, -1);
    labelListList ringEntities(4, labelList(m, -1));

    // Construct the hull
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
        hullEdges,
        hullFaces,
        hullCells,
        ringEntities
    );

    label numTriangulations = 0, isolatedVertex = -1;

    // Extract the appropriate triangulations
    extractTriangulation
    (
        0,
        m-1,
        K[checkIndex],
        numTriangulations,
        triangulations[checkIndex]
    );

    // Check old-volumes for the configuration.
    if
    (
        checkTriangulationVolumes
        (
            eIndex,
            hullVertices,
            triangulations[checkIndex]
        )
    )
    {
        // Reset all triangulations and bail out
        triangulations[checkIndex][0] = -1;
        triangulations[checkIndex][1] = -1;
        triangulations[checkIndex][2] = -1;

        return map;
    }

    // Determine the final swap triangulation
    label tF =
    (
        identify32Swap
        (
            eIndex,
            hullVertices,
            triangulations[checkIndex]
        )
    );

    // Check that the triangulation is valid

    if (tF == -1)
    {
        const edge& edgeToCheck = edges_[eIndex];

        Info << " All triangulations: " << nl
             << ' ' << triangulations[checkIndex][0] << nl
             << ' ' << triangulations[checkIndex][1] << nl
             << ' ' << triangulations[checkIndex][2] << nl
             << endl;

        FatalErrorIn
        (
            "\n"
            "const changeMap dynamicTopoFvMesh::removeEdgeFlips\n"
            "(\n"
            "    const label eIndex,\n"
            "    const scalar minQuality,\n"
            "    const PtrList<labelListList>& K,\n"
            "    PtrList<labelListList>& triangulations,\n"
            "    const label checkIndex\n"
            ")\n"
        )
            << "Could not determine 3-2 swap triangulation." << nl
            << "Edge: " << edgeToCheck << nl
            << "Edge Points: "
            << points_[edgeToCheck[0]] << ","
            << points_[edgeToCheck[1]] << nl
            << abort(FatalError);

        // Reset all triangulations and bail out
        triangulations[checkIndex][0] = -1;
        triangulations[checkIndex][1] = -1;
        triangulations[checkIndex][2] = -1;

        return map;
    }

    if (debug > 2)
    {
        Info << " Identified tF as: " << tF << endl;
        Info << " Triangulation: "
             << triangulations[checkIndex][0][tF] << " "
             << triangulations[checkIndex][1][tF] << " "
             << triangulations[checkIndex][2][tF] << " "
             << endl;
        Info << " All triangulations: " << nl
             << ' ' << triangulations[checkIndex][0] << nl
             << ' ' << triangulations[checkIndex][1] << nl
             << ' ' << triangulations[checkIndex][2] << nl
             << endl;
    }

    // Perform a series of 2-3 swaps
    label numSwaps = 0;

    while (numSwaps < (m-3))
    {
        for (label i = 0; i < (m-2); i++)
        {
            if ( (i != tF) && (triangulations[checkIndex][0][i] != -1) )
            {
                // Check if triangulation is on the boundary
                if
                (
                    boundaryTriangulation
                    (
                        i,
                        isolatedVertex,
                        triangulations[checkIndex]
                    )
                )
                {
                    // Perform 2-3 swap
                    map =
                    (
                        swap23
                        (
                            isolatedVertex,
                            eIndex,
                            i,
                            numTriangulations,
                            triangulations[checkIndex],
                            hullVertices,
                            hullFaces,
                            hullCells
                        )
                    );

                    // Done with this face, so reset it
                    triangulations[checkIndex][0][i] = -1;
                    triangulations[checkIndex][1][i] = -1;
                    triangulations[checkIndex][2][i] = -1;

                    numSwaps++;
                }
            }
        }

        if (numSwaps == 0)
        {
            Info << "Triangulations: " << endl;
            forAll(triangulations[checkIndex], row)
            {
                Info << triangulations[checkIndex][row] << endl;
            }

            // Should have performed at least one swap
            FatalErrorIn
            (
                "\n"
                "const changeMap dynamicTopoFvMesh::removeEdgeFlips\n"
                "(\n"
                "    const label eIndex,\n"
                "    const scalar minQuality,\n"
                "    const PtrList<labelListList>& K,\n"
                "    PtrList<labelListList>& triangulations,\n"
                "    const label checkIndex\n"
                ")\n"
            )
                << "Did not perform any 2-3 swaps" << nl
                << abort(FatalError);
        }
    }

    // Perform the final 3-2 / 2-2 swap
    map =
    (
        swap32
        (
            eIndex,
            tF,
            numTriangulations,
            triangulations[checkIndex],
            hullVertices,
            hullFaces,
            hullCells
        )
    );

    // Done with this face, so reset it
    triangulations[checkIndex][0][tF] = -1;
    triangulations[checkIndex][1][tF] = -1;
    triangulations[checkIndex][2][tF] = -1;

    // Finally remove the edge
    removeEdge(eIndex);

    // Increment the counter
    statistics_[1]++;

    // Set the flag
    topoChangeFlag_ = true;

    // Return a successful operation.
    return map;
}


// Extract triangulations from the programming table
void dynamicTopoFvMesh::extractTriangulation
(
    const label i,
    const label j,
    const labelListList& K,
    label& numTriangulations,
    labelListList& triangulations
) const
{
    if ( j >= (i+2) )
    {
        label k = K[i][j];

        // Fill in the triangulation list
        triangulations[0][numTriangulations] = i;
        triangulations[1][numTriangulations] = k;
        triangulations[2][numTriangulations] = j;

        // Increment triangulation count
        numTriangulations++;

        // Recursively call the function for the two sub-triangulations
        extractTriangulation(i,k,K,numTriangulations,triangulations);
        extractTriangulation(k,j,K,numTriangulations,triangulations);
    }
}


// Identify the 3-2 swap from the triangulation sequence
//  - Use an edge-plane intersection formula
label dynamicTopoFvMesh::identify32Swap
(
    const label eIndex,
    const labelList& hullVertices,
    const labelListList& triangulations
) const
{
    label m = hullVertices.size();
    const edge& edgeToCheck = edges_[eIndex];

    // Obtain intersection point.
    linePointRef segment
    (
        points_[edgeToCheck.start()],
        points_[edgeToCheck.end()]
    );

    vector intPt = vector::zero;

    // Configure a face with triangulation
    for (label i = 0; i < (m-2); i++)
    {
        bool intersects =
        (
            meshOps::segmentTriFaceIntersection
            (
                triPointRef
                (
                    points_[hullVertices[triangulations[0][i]]],
                    points_[hullVertices[triangulations[1][i]]],
                    points_[hullVertices[triangulations[2][i]]]
                ),
                segment,
                1e-6,
                intPt
            )
        );

        if (intersects)
        {
            return i;
        }
    }

    if (debug > 1)
    {
        Info << nl << nl << "Hull Vertices: " << endl;

        forAll(hullVertices, vertexI)
        {
            Info << hullVertices[vertexI] << ": "
                 << points_[hullVertices[vertexI]]
                 << endl;
        }

        InfoIn
        (
            "\n"
            "label dynamicTopoFvMesh::identify32Swap\n"
            "(\n"
            "    const label eIndex,\n"
            "    const labelList& hullVertices,\n"
            "    const labelListList& triangulations\n"
            ") const\n"
        )   << nl
            << "Could not determine 3-2 swap triangulation." << nl
            << "Edge: " << edgeToCheck << nl
            << "Edge Points: "
            << segment.start() << ","
            << segment.end() << nl
            << endl;
    }

    // Could not find an intersecting triangulation.
    //  - If this is a boundary edge, a curved surface-mesh.
    //    was probably the reason why this failed.
    //  - If so, declare the nearest triangulation instead.
    //  - Internal edges also occasionally encounter
    //    precision issues. Use the same approach.
    vector eCentre = segment.centre();

    scalarField dist(m-2, 0.0);

    label mT = -1, ePatch = whichEdgePatch(eIndex);
    bool foundTriangulation = false;

    for (label i = 0; i < (m-2); i++)
    {
        // Compute edge to face-centre distance.
        dist[i] =
        (
            mag
            (
                eCentre
              - triPointRef
                (
                    points_[hullVertices[triangulations[0][i]]],
                    points_[hullVertices[triangulations[1][i]]],
                    points_[hullVertices[triangulations[2][i]]]
                ).centre()
            )
        );
    }

    while (!foundTriangulation)
    {
        mT = findMin(dist);

        // Check validity for boundary edges
        if
        (
            (ePatch > -1) &&
            ((triangulations[0][mT] != 0) || (triangulations[2][mT] != m-1))
        )
        {
            // This is a 2-3 triangulation. Try again.
            dist[mT] = GREAT;
        }
        else
        {
            foundTriangulation = true;
        }
    }

    if (debug > 1)
    {
        Info << " All distances :" << dist << nl
             << " Triangulation index: " << mT
             << endl;
    }

    // Return the index
    return mT;
}


// Routine to check whether the triangulation at the
// index lies on the boundary of the vertex ring.
bool dynamicTopoFvMesh::boundaryTriangulation
(
    const label index,
    label& isolatedVertex,
    labelListList& triangulations
) const
{
    label first = 0, second = 0, third = 0;

    // Count for occurrences
    forAll(triangulations, row)
    {
        forAll(triangulations[row], col)
        {
            if (triangulations[row][col] == triangulations[0][index])
            {
                first++;
            }

            if (triangulations[row][col] == triangulations[1][index])
            {
                second++;
            }

            if (triangulations[row][col] == triangulations[2][index])
            {
                third++;
            }
        }
    }

    if (first == 1)
    {
        isolatedVertex = triangulations[0][index];
        return true;
    }

    if (second == 1)
    {
        isolatedVertex = triangulations[1][index];
        return true;
    }

    if (third == 1)
    {
        isolatedVertex = triangulations[2][index];
        return true;
    }

    // This isn't a boundary triangulation
    return false;
}


// Utility method to compute the minimum quality of a vertex hull
scalar dynamicTopoFvMesh::computeMinQuality
(
    const label eIndex
) const
{
    scalar minQuality = GREAT;
    scalar cQuality = 0.0;

    // Obtain a reference to this edge and corresponding edgePoints
    const edge& edgeToCheck = edges_[eIndex];
    const labelList& hullVertices = edgePoints_[eIndex];

    // Obtain point references
    const point& a = points_[edgeToCheck[0]];
    const point& c = points_[edgeToCheck[1]];

    if (whichEdgePatch(eIndex) < 0)
    {
        // Internal edge.
        forAll(hullVertices, indexI)
        {
            label prevIndex = hullVertices.rcIndex(indexI);

            // Pick vertices off the list
            const point& b = points_[hullVertices[prevIndex]];
            const point& d = points_[hullVertices[indexI]];

            // Compute the quality
            cQuality = tetMetric_(a, b, c, d);

            // Check if the quality is worse
            minQuality = Foam::min(cQuality, minQuality);
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

            // Compute the quality
            cQuality = tetMetric_(a, b, c, d);

            // Check if the quality is worse
            minQuality = Foam::min(cQuality, minQuality);
        }
    }

    // Ensure that the mesh is valid
    if (minQuality < 0.0)
    {
        // if (debug > 3)
        {
            // Write out faces and cells for post processing.
            labelHashSet iFaces, iCells, bFaces;

            const labelList& eFaces = edgeFaces_[eIndex];

            forAll(eFaces, faceI)
            {
                iFaces.insert(eFaces[faceI]);

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
            writeVTK(Foam::name(eIndex) + "_iFaces", iFaces.toc(), 2);

            // Write out the boundary patches (for post-processing reference)
            for
            (
                label faceI = nOldInternalFaces_;
                faceI < faces_.size();
                faceI++
            )
            {
                if (faces_[faceI].empty())
                {
                    continue;
                }

                label pIndex = whichPatch(faceI);

                if (pIndex != -1)
                {
                    bFaces.insert(faceI);
                }
            }

            writeVTK(Foam::name(eIndex) + "_bFaces", bFaces.toc(), 2);
        }

        FatalErrorIn
        (
            "scalar dynamicTopoFvMesh::computeMinQuality"
            "(const label eIndex) const"
        )
            << "Encountered negative cell-quality!" << nl
            << "Edge: " << eIndex << ": " << edgeToCheck << nl
            << "EdgePoints: " << hullVertices << nl
            << "Minimum Quality: " << minQuality
            << abort(FatalError);
    }

    return minQuality;
}


// Method used to perform a 2-3 swap in 3D
// - Returns a changeMap with the index of
//   the triangulated face in opposingFace.
const changeMap dynamicTopoFvMesh::swap23
(
    const label isolatedVertex,
    const label eIndex,
    const label triangulationIndex,
    const label numTriangulations,
    const labelListList& triangulations,
    const labelList& hullVertices,
    const labelList& hullFaces,
    const labelList& hullCells
)
{
    // A 2-3 swap performs the following operations:
    //      [1] Remove face: [ edge[0] edge[1] isolatedVertex ]
    //      [2] Remove two cells on either side of removed face
    //      [3] Add one edge
    //      [4] Add three new faces
    //      [5] Add three new cells
    //      Update faceEdges, edgeFaces and edgePoints information

    changeMap map;

    // Obtain a copy of the edge
    edge edgeToCheck = edges_[eIndex];

    label faceForRemoval = hullFaces[isolatedVertex];
    label vertexForRemoval = hullVertices[isolatedVertex];

    // Determine the two cells to be removed
    FixedList<label,2> cellsForRemoval;
    cellsForRemoval[0] = owner_[faceForRemoval];
    cellsForRemoval[1] = neighbour_[faceForRemoval];

    if (debug > 1)
    {
        // Print out arguments
        Info << endl;
        Info << "== Swapping 2-3 ==" << endl;
        Info << "Edge: " << eIndex << ": " << edgeToCheck << endl;

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

            Info << "Ring: " << hullVertices << endl;
            Info << "Faces: " << hullFaces << endl;
            Info << "Cells: " << hullCells << endl;
            Info << "Triangulation: "
                 << triangulations[0][triangulationIndex] << " "
                 << triangulations[1][triangulationIndex] << " "
                 << triangulations[2][triangulationIndex] << " "
                 << endl;
            Info << "Isolated vertex: " << isolatedVertex << endl;
        }

        if (debug > 3)
        {
            writeVTK
            (
                Foam::name(eIndex)
              + '(' + Foam::name(edgeToCheck[0])
              + ',' + Foam::name(edgeToCheck[1]) + ')'
              + "_beforeSwap_"
              + Foam::name(numTriangulations) + "_"
              + Foam::name(triangulationIndex),
                cellsForRemoval
            );
        }
    }

    // Check if this is an internal face
    if (cellsForRemoval[1] == -1)
    {
        // Write out for post-processing
        writeVTK("Edge23_" + Foam::name(eIndex), eIndex, 1);
        writeVTK("Cells23_" + Foam::name(eIndex), hullCells, 3);

        FatalErrorIn
        (
            "\n"
            "const changeMap dynamicTopoFvMesh::swap23\n"
            "(\n"
            "    const label isolatedVertex,\n"
            "    const label eIndex,\n"
            "    const label triangulationIndex,\n"
            "    const label numTriangulations,\n"
            "    const labelListList& triangulations,\n"
            "    const labelList& hullVertices,\n"
            "    const labelList& hullFaces,\n"
            "    const labelList& hullCells\n"
            ")\n"
        )
            << " Expected an internal face,"
            << " but found a boundary one instead." << nl
            << " Looks like identify32Swap couldn't correctly identify"
            << " the 2-2 swap triangulation." << nl
            << abort(FatalError);
    }

    // Add three new cells to the end of the cell list
    FixedList<label,3> newCellIndex(-1);
    FixedList<cell, 3> newTetCell(cell(4));

    forAll(newCellIndex, cellI)
    {
        scalar avgScale = -1.0;

        if (edgeRefinement_)
        {
            avgScale =
            (
                0.5 *
                (
                    lengthScale_[cellsForRemoval[0]]
                  + lengthScale_[cellsForRemoval[1]]
                )
            );
        }

        // Insert the cell
        newCellIndex[cellI] = insertCell(newTetCell[cellI], avgScale);
    }

    // Obtain point-ordering for the other vertices
    // otherVertices[0] is the point before isolatedVertex
    // otherVertices[1] is the point after isolatedVertex
    FixedList<label,2> otherVertices;

    if (triangulations[0][triangulationIndex] == isolatedVertex)
    {
        otherVertices[0] = hullVertices[triangulations[2][triangulationIndex]];
        otherVertices[1] = hullVertices[triangulations[1][triangulationIndex]];
    }
    else
    if (triangulations[1][triangulationIndex] == isolatedVertex)
    {
        otherVertices[0] = hullVertices[triangulations[0][triangulationIndex]];
        otherVertices[1] = hullVertices[triangulations[2][triangulationIndex]];
    }
    else
    {
        otherVertices[0] = hullVertices[triangulations[1][triangulationIndex]];
        otherVertices[1] = hullVertices[triangulations[0][triangulationIndex]];
    }

    // Insert three new internal faces
    FixedList<label,3> newFaceIndex;
    face tmpTriFace(3);

    // First face: The actual triangulation
    tmpTriFace[0] = otherVertices[0];
    tmpTriFace[1] = vertexForRemoval;
    tmpTriFace[2] = otherVertices[1];

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

    // Add an entry to the map
    map.opposingFace() = newFaceIndex[0];

    // Second face: Triangle involving edgeToCheck[0]
    tmpTriFace[0] = otherVertices[0];
    tmpTriFace[1] = edgeToCheck[0];
    tmpTriFace[2] = otherVertices[1];

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

    // Third face: Triangle involving edgeToCheck[1]
    tmpTriFace[0] = otherVertices[1];
    tmpTriFace[1] = edgeToCheck[1];
    tmpTriFace[2] = otherVertices[0];

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

    // Append three dummy faceEdges entries.
    for (label i = 0; i < 3; i++)
    {
        faceEdges_.append(labelList(0));
    }

    // Add an entry to edgeFaces
    labelList newEdgeFaces(3);
    newEdgeFaces[0] = newFaceIndex[0];
    newEdgeFaces[1] = newFaceIndex[1];
    newEdgeFaces[2] = newFaceIndex[2];

    // Add an entry for edgePoints as well
    labelList newEdgePoints(3);
    newEdgePoints[0] = vertexForRemoval;
    newEdgePoints[1] = edgeToCheck[0];
    newEdgePoints[2] = edgeToCheck[1];

    // Add a new internal edge to the mesh
    label newEdgeIndex =
    (
        insertEdge
        (
            -1,
            edge
            (
                otherVertices[0],
                otherVertices[1]
            ),
            newEdgeFaces,
            newEdgePoints
        )
    );

    // Define the six edges to check while building faceEdges:
    FixedList<edge,6> check;

    check[0][0] = vertexForRemoval; check[0][1] = otherVertices[0];
    check[1][0] = vertexForRemoval; check[1][1] = otherVertices[1];
    check[2][0] = edgeToCheck[0];   check[2][1] = otherVertices[0];
    check[3][0] = edgeToCheck[1];   check[3][1] = otherVertices[0];
    check[4][0] = edgeToCheck[0];   check[4][1] = otherVertices[1];
    check[5][0] = edgeToCheck[1];   check[5][1] = otherVertices[1];

    // Add three new entries to faceEdges
    label nE0 = 0, nE1 = 0, nE2 = 0;
    FixedList<labelList,3> newFaceEdges(labelList(3));

    newFaceEdges[0][nE0++] = newEdgeIndex;
    newFaceEdges[1][nE1++] = newEdgeIndex;
    newFaceEdges[2][nE2++] = newEdgeIndex;

    // Fill-in information for the three new cells,
    // and correct info on existing neighbouring cells
    label nF0 = 0, nF1 = 0, nF2 = 0;
    FixedList<bool,2> foundEdge;

    // Add the newly created faces to cells
    newTetCell[0][nF0++] = newFaceIndex[0];
    newTetCell[0][nF0++] = newFaceIndex[2];
    newTetCell[1][nF1++] = newFaceIndex[0];
    newTetCell[1][nF1++] = newFaceIndex[1];
    newTetCell[2][nF2++] = newFaceIndex[1];
    newTetCell[2][nF2++] = newFaceIndex[2];

    forAll(cellsForRemoval, cellI)
    {
        label cellIndex = cellsForRemoval[cellI];

        forAll(cells_[cellIndex], faceI)
        {
            label faceIndex = cells_[cellIndex][faceI];

            foundEdge[0] = false; foundEdge[1] = false;

            // Check if face contains edgeToCheck[0]
            if
            (
                (faces_[faceIndex][0] == edgeToCheck[0])
             || (faces_[faceIndex][1] == edgeToCheck[0])
             || (faces_[faceIndex][2] == edgeToCheck[0])
            )
            {
                foundEdge[0] = true;
            }

            // Check if face contains edgeToCheck[1]
            if
            (
                (faces_[faceIndex][0] == edgeToCheck[1])
             || (faces_[faceIndex][1] == edgeToCheck[1])
             || (faces_[faceIndex][2] == edgeToCheck[1])
            )
            {
                foundEdge[1] = true;
            }

            // Face is connected to edgeToCheck[0]
            if (foundEdge[0] && !foundEdge[1])
            {
                // Check if a face-flip is necessary
                if (owner_[faceIndex] == cellIndex)
                {
                    if (neighbour_[faceIndex] == -1)
                    {
                        // Change the owner
                        owner_[faceIndex] = newCellIndex[1];
                    }
                    else
                    {
                        // Flip this face
                        faces_[faceIndex] = faces_[faceIndex].reverseFace();
                        owner_[faceIndex] = neighbour_[faceIndex];
                        neighbour_[faceIndex] = newCellIndex[1];

                        setFlip(faceIndex);
                    }
                }
                else
                {
                    // Flip is unnecessary. Just update neighbour
                    neighbour_[faceIndex] = newCellIndex[1];
                }

                // Add this face to the cell
                newTetCell[1][nF1++] = faceIndex;

                // Update faceEdges, edgeFaces, and edgePoints.
                forAll(faceEdges_[faceIndex], edgeI)
                {
                    if (edges_[faceEdges_[faceIndex][edgeI]] == check[0])
                    {
                        newFaceEdges[0][nE0++] = faceEdges_[faceIndex][edgeI];

                        meshOps::sizeUpList
                        (
                            newFaceIndex[0],
                            edgeFaces_[faceEdges_[faceIndex][edgeI]]
                        );

                        meshOps::insertLabel
                        (
                            otherVertices[1],
                            edgeToCheck[0],
                            edgeToCheck[1],
                            edgePoints_[faceEdges_[faceIndex][edgeI]]
                        );
                    }

                    if (edges_[faceEdges_[faceIndex][edgeI]] == check[1])
                    {
                        newFaceEdges[0][nE0++] = faceEdges_[faceIndex][edgeI];

                        meshOps::sizeUpList
                        (
                            newFaceIndex[0],
                            edgeFaces_[faceEdges_[faceIndex][edgeI]]
                        );

                        meshOps::insertLabel
                        (
                            otherVertices[0],
                            edgeToCheck[0],
                            edgeToCheck[1],
                            edgePoints_[faceEdges_[faceIndex][edgeI]]
                        );
                    }

                    if (edges_[faceEdges_[faceIndex][edgeI]] == check[2])
                    {
                        newFaceEdges[1][nE1++] = faceEdges_[faceIndex][edgeI];

                        meshOps::sizeUpList
                        (
                            newFaceIndex[1],
                            edgeFaces_[faceEdges_[faceIndex][edgeI]]
                        );

                        meshOps::insertLabel
                        (
                            otherVertices[1],
                            vertexForRemoval,
                            edgeToCheck[1],
                            edgePoints_[faceEdges_[faceIndex][edgeI]]
                        );
                    }

                    if (edges_[faceEdges_[faceIndex][edgeI]] == check[4])
                    {
                        newFaceEdges[1][nE1++] = faceEdges_[faceIndex][edgeI];

                        meshOps::sizeUpList
                        (
                            newFaceIndex[1],
                            edgeFaces_[faceEdges_[faceIndex][edgeI]]
                        );

                        meshOps::insertLabel
                        (
                            otherVertices[0],
                            vertexForRemoval,
                            edgeToCheck[1],
                            edgePoints_[faceEdges_[faceIndex][edgeI]]
                        );
                    }
                }
            }

            // Face is connected to edgeToCheck[1]
            if (!foundEdge[0] && foundEdge[1])
            {
                // Check if a face-flip is necessary
                if (owner_[faceIndex] == cellIndex)
                {
                    if (neighbour_[faceIndex] == -1)
                    {
                        // Change the owner
                        owner_[faceIndex] = newCellIndex[0];
                    }
                    else
                    {
                        // Flip this face
                        faces_[faceIndex] = faces_[faceIndex].reverseFace();
                        owner_[faceIndex] = neighbour_[faceIndex];
                        neighbour_[faceIndex] = newCellIndex[0];

                        setFlip(faceIndex);
                    }
                }
                else
                {
                    // Flip is unnecessary. Just update neighbour
                    neighbour_[faceIndex] = newCellIndex[0];
                }

                // Add this face to the cell
                newTetCell[0][nF0++] = faceIndex;

                // Update faceEdges, edgeFaces, and edgePoints.
                const labelList& fEdges = faceEdges_[faceIndex];

                forAll(fEdges, edgeI)
                {
                    if (edges_[fEdges[edgeI]] == check[3])
                    {
                        newFaceEdges[2][nE2++] = fEdges[edgeI];

                        meshOps::sizeUpList
                        (
                            newFaceIndex[2],
                            edgeFaces_[fEdges[edgeI]]
                        );

                        meshOps::insertLabel
                        (
                            otherVertices[1],
                            vertexForRemoval,
                            edgeToCheck[0],
                            edgePoints_[fEdges[edgeI]]
                        );
                    }

                    if (edges_[fEdges[edgeI]] == check[5])
                    {
                        newFaceEdges[2][nE2++] = fEdges[edgeI];

                        meshOps::sizeUpList
                        (
                            newFaceIndex[2],
                            edgeFaces_[fEdges[edgeI]]
                        );

                        meshOps::insertLabel
                        (
                            otherVertices[0],
                            vertexForRemoval,
                            edgeToCheck[0],
                            edgePoints_[fEdges[edgeI]]
                        );
                    }
                }
            }

            // Face is connected to both edgeToCheck [0] and [1]
            if
            (
                (foundEdge[0] && foundEdge[1]) &&
                (faceIndex != faceForRemoval)
            )
            {
                // Check if a face-flip is necessary
                if (owner_[faceIndex] == cellIndex)
                {
                    if (neighbour_[faceIndex] == -1)
                    {
                        // Change the owner
                        owner_[faceIndex] = newCellIndex[2];
                    }
                    else
                    {
                        // Flip this face
                        faces_[faceIndex] = faces_[faceIndex].reverseFace();
                        owner_[faceIndex] = neighbour_[faceIndex];
                        neighbour_[faceIndex] = newCellIndex[2];

                        setFlip(faceIndex);
                    }
                }
                else
                {
                    // Flip is unnecessary. Just update neighbour
                    neighbour_[faceIndex] = newCellIndex[2];
                }

                // Add this face to the cell
                newTetCell[2][nF2++] = faceIndex;
            }
        }
    }

    // Now update faceEdges for the three new faces
    forAll(newFaceEdges, faceI)
    {
        faceEdges_[newFaceIndex[faceI]] = newFaceEdges[faceI];
    }

    // Update edgeFaces and edgePoints for edges of the removed face
    label otherPoint = -1, nextPoint = -1;

    forAll(faceEdges_[faceForRemoval], edgeI)
    {
        label edgeIndex = faceEdges_[faceForRemoval][edgeI];

        meshOps::sizeDownList
        (
            faceForRemoval,
            edgeFaces_[edgeIndex]
        );

        // Find the isolated point and remove it
        meshOps::findIsolatedPoint
        (
            faces_[faceForRemoval],
            edges_[edgeIndex],
            otherPoint,
            nextPoint
        );

        meshOps::sizeDownList
        (
            otherPoint,
            edgePoints_[edgeIndex]
        );
    }

    // Remove the face
    removeFace(faceForRemoval);

    forAll(cellsForRemoval, cellI)
    {
        removeCell(cellsForRemoval[cellI]);
    }

    // Fill-in candidate mapping information
    labelList mC(2, -1);

    forAll(mC, indexI)
    {
        mC[indexI] = cellsForRemoval[indexI];
    }

    // Update the cell list with newly configured cells.
    forAll(newCellIndex, cellI)
    {
        cells_[newCellIndex[cellI]] = newTetCell[cellI];

        if (cellI == 2)
        {
            // Skip mapping for the intermediate cell.
            setCellMapping(newCellIndex[cellI], mC, false);
        }
        else
        {
            // Set the mapping for this cell
            setCellMapping(newCellIndex[cellI], mC);
        }
    }

    // Fill in mapping information for three new faces.
    // Since they're all internal, interpolate fluxes by default.
    forAll(newFaceIndex, faceI)
    {
        setFaceMapping(newFaceIndex[faceI]);
    }

    if (debug > 2)
    {
        Info << "Added edge: " << endl;

        Info << newEdgeIndex << ":: "
             << edges_[newEdgeIndex]
             << " edgeFaces: "
             << edgeFaces_[newEdgeIndex]
             << endl;

        Info << "Added faces: " << endl;

        forAll(newFaceIndex, faceI)
        {
            Info << newFaceIndex[faceI] << ":: "
                 << faces_[newFaceIndex[faceI]]
                 << " faceEdges: "
                 << faceEdges_[newFaceIndex[faceI]]
                 << endl;
        }

        Info << "Added cells: " << endl;

        forAll(newCellIndex, cellI)
        {
            Info << newCellIndex[cellI] << ":: "
                 << cells_[newCellIndex[cellI]]
                 << endl;
        }

        if (debug > 3)
        {
            writeVTK
            (
                Foam::name(eIndex)
              + '(' + Foam::name(edgeToCheck[0])
              + ',' + Foam::name(edgeToCheck[1]) + ')'
              + "_afterSwap_"
              + Foam::name(numTriangulations) + "_"
              + Foam::name(triangulationIndex),
                newCellIndex
            );
        }
    }

    // Specify that the swap was successful
    map.type() = 1;

    // Return the changeMap
    return map;
}


// Method used to perform a 2-2 / 3-2 swap in 3D
// - Returns a changeMap with a type specifying:
//     1: Swap was successful
// - The index of the triangulated face in opposingFace.
// - For 2-2 swaps, the newly added edge index.
const changeMap dynamicTopoFvMesh::swap32
(
    const label eIndex,
    const label triangulationIndex,
    const label numTriangulations,
    const labelListList& triangulations,
    const labelList& hullVertices,
    const labelList& hullFaces,
    const labelList& hullCells
)
{
    // A 2-2 / 3-2 swap performs the following operations:
    //      [1] Remove three faces surrounding edgeToCheck
    //      [2] Remove two (2-2 swap) or three(3-2 swap)
    //          cells surrounding edgeToCheck
    //      [3] Add one internal face
    //      [4] Add two new cells
    //      [5] If edgeToCheck is on a boundary,
    //          add two boundary faces and a boundary edge (2-2 swap)
    //      eIndex is removed later by removeEdgeFlips
    //      Update faceEdges, edgeFaces and edgePoints information

    changeMap map;

    // Obtain a copy of the edge
    edge edgeToCheck = edges_[eIndex];

    // Determine the patch this edge belongs to
    label edgePatch = whichEdgePatch(eIndex);

    // Determine the three faces to be removed
    FixedList<label,3> facesForRemoval;
    labelHashSet cellsForRemoval(3);

    forAll(facesForRemoval, faceI)
    {
        facesForRemoval[faceI] =
        (
            hullFaces[triangulations[faceI][triangulationIndex]]
        );

        label own = owner_[facesForRemoval[faceI]];
        label nei = neighbour_[facesForRemoval[faceI]];

        // Check and add cells as well
        if (!cellsForRemoval.found(own))
        {
            cellsForRemoval.insert(own);
        }

        if (!cellsForRemoval.found(nei) && nei != -1)
        {
            cellsForRemoval.insert(nei);
        }
    }

    labelList cellRemovalList = cellsForRemoval.toc();

    if (debug > 1)
    {
        // Print out arguments
        Info << endl;

        if (edgePatch < 0)
        {
            Info << "== Swapping 3-2 ==" << endl;
        }
        else
        {
            Info << "== Swapping 2-2 ==" << endl;
        }

        Info << "Edge: " << eIndex << ": " << edgeToCheck << endl;

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

            Info << "Ring: " << hullVertices << endl;
            Info << "Faces: " << hullFaces << endl;
            Info << "Cells: " << hullCells << endl;
            Info << "Triangulation: "
                 << triangulations[0][triangulationIndex] << " "
                 << triangulations[1][triangulationIndex] << " "
                 << triangulations[2][triangulationIndex] << " "
                 << endl;
        }

        if (debug > 3)
        {
            writeVTK
            (
                Foam::name(eIndex)
              + '(' + Foam::name(edgeToCheck[0])
              + ',' + Foam::name(edgeToCheck[1]) + ')'
              + "_beforeSwap_"
              + Foam::name(numTriangulations) + "_"
              + Foam::name(triangulationIndex),
                cellRemovalList
            );
        }
    }

    // Add two new cells to the end of the cell list
    FixedList<label,2> newCellIndex(-1);
    FixedList<cell, 2> newTetCell(cell(4));

    forAll(newCellIndex, cellI)
    {
        scalar avgScale = 0.0;

        if (edgeRefinement_)
        {
            forAll(cellRemovalList, indexI)
            {
                avgScale += lengthScale_[cellRemovalList[indexI]];
            }

            avgScale /= cellRemovalList.size();
        }

        // Insert the cell
        newCellIndex[cellI] = insertCell(newTetCell[cellI], avgScale);
    }

    // Insert a new internal face
    face newTriFace(3);

    newTriFace[0] = hullVertices[triangulations[0][triangulationIndex]];
    newTriFace[1] = hullVertices[triangulations[1][triangulationIndex]];
    newTriFace[2] = hullVertices[triangulations[2][triangulationIndex]];

    label newFaceIndex =
    (
        insertFace
        (
            -1,
            newTriFace,
            newCellIndex[0],
            newCellIndex[1]
        )
    );

    // Add an entry to the map
    map.opposingFace() = newFaceIndex;

    // Add faceEdges for the new face as well.
    faceEdges_.append(labelList(3));

    // Define the three edges to check while building faceEdges:
    FixedList<edge,3> check;

    check[0][0] = newTriFace[0]; check[0][1] = newTriFace[1];
    check[1][0] = newTriFace[1]; check[1][1] = newTriFace[2];
    check[2][0] = newTriFace[2]; check[2][1] = newTriFace[0];

    // For 2-2 swaps, two faces are introduced
    label nE = 0, nBf = 0;
    FixedList<label,2> nBE(0);
    FixedList<labelList,2> bdyFaceEdges(labelList(3, -1));

    // Fill-in information for the two new cells,
    // and correct info on existing neighbouring cells
    label nF0 = 0, nF1 = 0;
    label otherPoint = -1, nextPoint = -1;
    FixedList<bool,2> foundEdge;

    // For a 2-2 swap on a boundary edge,
    // add two boundary faces and an edge
    label newEdgeIndex = -1;
    labelList oldBdyFaceIndex(2, -1), newBdyFaceIndex(2, -1);

    if (edgePatch > -1)
    {
        // Temporary local variables
        label facePatch = -1;
        edge newEdge(-1, -1);
        FixedList<label,2> nBEdge(0);
        FixedList<FixedList<label,2>,2> bdyEdges;
        FixedList<face,2> newBdyTriFace(face(3));

        // Get a cue for face orientation from existing faces
        forAll(facesForRemoval, faceI)
        {
            if (neighbour_[facesForRemoval[faceI]] == -1)
            {
                facePatch = whichPatch(facesForRemoval[faceI]);

                // Record this face-index for mapping.
                oldBdyFaceIndex[nBf++] = facesForRemoval[faceI];

                meshOps::findIsolatedPoint
                (
                    faces_[facesForRemoval[faceI]],
                    edgeToCheck,
                    otherPoint,
                    nextPoint
                );

                if (nextPoint == edgeToCheck[0])
                {
                    newEdge[1] = otherPoint;
                    newBdyTriFace[0][0] = otherPoint;
                    newBdyTriFace[0][1] = edgeToCheck[0];
                    newBdyTriFace[1][2] = otherPoint;
                }
                else
                {
                    newEdge[0] = otherPoint;
                    newBdyTriFace[1][0] = otherPoint;
                    newBdyTriFace[1][1] = edgeToCheck[1];
                    newBdyTriFace[0][2] = otherPoint;
                }

                // Also update faceEdges for the new boundary faces
                forAll(faceEdges_[facesForRemoval[faceI]], edgeI)
                {
                    if
                    (
                        edges_[faceEdges_[facesForRemoval[faceI]][edgeI]]
                     == edge(edgeToCheck[0], otherPoint)
                    )
                    {
                        bdyFaceEdges[0][nBE[0]++] =
                        (
                            faceEdges_[facesForRemoval[faceI]][edgeI]
                        );

                        bdyEdges[0][nBEdge[0]++] =
                        (
                            faceEdges_[facesForRemoval[faceI]][edgeI]
                        );
                    }

                    if
                    (
                        edges_[faceEdges_[facesForRemoval[faceI]][edgeI]]
                     == edge(edgeToCheck[1], otherPoint)
                    )
                    {
                        bdyFaceEdges[1][nBE[1]++] =
                        (
                            faceEdges_[facesForRemoval[faceI]][edgeI]
                        );

                        bdyEdges[1][nBEdge[1]++] =
                        (
                            faceEdges_[facesForRemoval[faceI]][edgeI]
                        );
                    }
                }
            }
        }

        // Insert the first of two new faces
        newBdyFaceIndex[0] =
        (
            insertFace
            (
                facePatch,
                newBdyTriFace[0],
                newCellIndex[1],
                -1
            )
        );

        // Add this face to the map.
        map.addFace(newBdyFaceIndex[0]);

        // Insert the second of two new faces
        newBdyFaceIndex[1] =
        (
            insertFace
            (
                facePatch,
                newBdyTriFace[1],
                newCellIndex[0],
                -1
            )
        );

        // Add this face to the map.
        map.addFace(newBdyFaceIndex[1]);

        // Update the new cells
        newTetCell[0][nF0++] = newBdyFaceIndex[1];
        newTetCell[1][nF1++] = newBdyFaceIndex[0];

        // Add an edgeFaces entry
        labelList newBdyEdgeFaces(3, -1);
        newBdyEdgeFaces[0] = newBdyFaceIndex[0];
        newBdyEdgeFaces[1] = newFaceIndex;
        newBdyEdgeFaces[2] = newBdyFaceIndex[1];

        // Find the point other than the new edge
        // on the new triangular face
        meshOps::findIsolatedPoint
        (
            newTriFace,
            newEdge,
            otherPoint,
            nextPoint
        );

        // Add an edgePoints entry
        labelList newBdyEdgePoints(3, -1);
        newBdyEdgePoints[0] = edgeToCheck[0];
        newBdyEdgePoints[1] = otherPoint;
        newBdyEdgePoints[2] = edgeToCheck[1];

        // Insert the edge
        newEdgeIndex =
        (
            insertEdge
            (
                edgePatch,
                newEdge,
                newBdyEdgeFaces,
                newBdyEdgePoints
            )
        );

        // Add this edge to the map.
        map.addEdge(newEdgeIndex);

        // Update faceEdges with the new edge
        faceEdges_[newFaceIndex][nE++] = newEdgeIndex;
        bdyFaceEdges[0][nBE[0]++] = newEdgeIndex;
        bdyFaceEdges[1][nBE[1]++] = newEdgeIndex;

        // Update edgeFaces and edgePoints with the two new faces
        forAll(bdyEdges[0], edgeI)
        {
            meshOps::sizeUpList
            (
                newBdyFaceIndex[0],
                edgeFaces_[bdyEdges[0][edgeI]]
            );

            meshOps::sizeUpList
            (
                newBdyFaceIndex[1],
                edgeFaces_[bdyEdges[1][edgeI]]
            );

            // Replace the edgePoints label, and preserve position on the list
            meshOps::findIsolatedPoint
            (
                newBdyTriFace[0],
                edges_[bdyEdges[0][edgeI]],
                otherPoint,
                nextPoint
            );

            meshOps::replaceLabel
            (
                edgeToCheck[1],
                otherPoint,
                edgePoints_[bdyEdges[0][edgeI]]
            );

            // Size up edgePoints again, so that it is sized down later
            meshOps::sizeUpList
            (
                edgeToCheck[1],
                edgePoints_[bdyEdges[0][edgeI]]
            );

            // Replace the edgePoints label, and preserve position on the list
            meshOps::findIsolatedPoint
            (
                newBdyTriFace[1],
                edges_[bdyEdges[1][edgeI]],
                otherPoint,
                nextPoint
            );

            meshOps::replaceLabel
            (
                edgeToCheck[0],
                otherPoint,
                edgePoints_[bdyEdges[1][edgeI]]
            );

            // Size up edgePoints again, so that it is sized down later
            meshOps::sizeUpList
            (
                edgeToCheck[0],
                edgePoints_[bdyEdges[1][edgeI]]
            );
        }

        // Add faceEdges for the two new boundary faces
        faceEdges_.append(bdyFaceEdges[0]);
        faceEdges_.append(bdyFaceEdges[1]);

        // Update the number of surface swaps.
        statistics_[2]++;
    }

    newTetCell[0][nF0++] = newFaceIndex;
    newTetCell[1][nF1++] = newFaceIndex;

    forAll(cellRemovalList, cellI)
    {
        label cellIndex = cellRemovalList[cellI];

        forAll(cells_[cellIndex], faceI)
        {
            label faceIndex = cells_[cellIndex][faceI];

            foundEdge[0] = false; foundEdge[1] = false;

            // Find the face that contains only
            // edgeToCheck[0] or edgeToCheck[1]
            forAll(faces_[faceIndex], pointI)
            {
                if (faces_[faceIndex][pointI] == edgeToCheck[0])
                {
                    foundEdge[0] = true;
                }

                if (faces_[faceIndex][pointI] == edgeToCheck[1])
                {
                    foundEdge[1] = true;
                }
            }

            // Face is connected to edgeToCheck[0]
            if (foundEdge[0] && !foundEdge[1])
            {
                // Check if a face-flip is necessary
                if (owner_[faceIndex] == cellIndex)
                {
                    if (neighbour_[faceIndex] == -1)
                    {
                        // Change the owner
                        owner_[faceIndex] = newCellIndex[1];
                    }
                    else
                    {
                        // Flip this face
                        faces_[faceIndex] = faces_[faceIndex].reverseFace();
                        owner_[faceIndex] = neighbour_[faceIndex];
                        neighbour_[faceIndex] = newCellIndex[1];

                        setFlip(faceIndex);
                    }
                }
                else
                {
                    // Flip is unnecessary. Just update neighbour
                    neighbour_[faceIndex] = newCellIndex[1];
                }

                // Add this face to the cell
                newTetCell[1][nF1++] = faceIndex;

                // Update faceEdges, edgeFaces and edgePoints
                forAll(faceEdges_[faceIndex], edgeI)
                {
                    if
                    (
                        (edges_[faceEdges_[faceIndex][edgeI]] == check[0])
                     || (edges_[faceEdges_[faceIndex][edgeI]] == check[1])
                     || (edges_[faceEdges_[faceIndex][edgeI]] == check[2])
                    )
                    {
                        faceEdges_[newFaceIndex][nE++] =
                        (
                            faceEdges_[faceIndex][edgeI]
                        );

                        meshOps::sizeUpList
                        (
                            newFaceIndex,
                            edgeFaces_[faceEdges_[faceIndex][edgeI]]
                        );

                        // Find the isolated point and insert it
                        meshOps::findIsolatedPoint
                        (
                            newTriFace,
                            edges_[faceEdges_[faceIndex][edgeI]],
                            otherPoint,
                            nextPoint
                        );

                        meshOps::insertLabel
                        (
                            otherPoint,
                            edgeToCheck[0],
                            edgeToCheck[1],
                            edgePoints_[faceEdges_[faceIndex][edgeI]]
                        );

                        break;
                    }
                }
            }

            // Face is connected to edgeToCheck[1]
            if (!foundEdge[0] && foundEdge[1])
            {
                // Check if a face-flip is necessary
                if (owner_[faceIndex] == cellIndex)
                {
                    if (neighbour_[faceIndex] == -1)
                    {
                        // Change the owner
                        owner_[faceIndex] = newCellIndex[0];
                    }
                    else
                    {
                        // Flip this face
                        faces_[faceIndex] = faces_[faceIndex].reverseFace();
                        owner_[faceIndex] = neighbour_[faceIndex];
                        neighbour_[faceIndex] = newCellIndex[0];

                        setFlip(faceIndex);
                    }
                }
                else
                {
                    // Flip is unnecessary. Just update neighbour
                    neighbour_[faceIndex] = newCellIndex[0];
                }

                // Add this face to the cell
                newTetCell[0][nF0++] = faceIndex;
            }
        }
    }

    // Remove the faces and update associated edges
    forAll(facesForRemoval, faceI)
    {
        // Update edgeFaces and edgePoints
        forAll(faceEdges_[facesForRemoval[faceI]], edgeI)
        {
            label edgeIndex = faceEdges_[facesForRemoval[faceI]][edgeI];

            if (edgeIndex != eIndex)
            {
                meshOps::sizeDownList
                (
                    facesForRemoval[faceI],
                    edgeFaces_[edgeIndex]
                );

                // Find the isolated point and remove it
                meshOps::findIsolatedPoint
                (
                    faces_[facesForRemoval[faceI]],
                    edges_[edgeIndex],
                    otherPoint,
                    nextPoint
                );

                meshOps::sizeDownList
                (
                    otherPoint,
                    edgePoints_[edgeIndex]
                );
            }
        }

        // Now remove the face...
        removeFace(facesForRemoval[faceI]);
    }

    forAll(cellRemovalList, cellI)
    {
        removeCell(cellRemovalList[cellI]);
    }

    // Update the cell list with newly configured cells.
    forAll(newCellIndex, cellI)
    {
        cells_[newCellIndex[cellI]] = newTetCell[cellI];

        // Fill-in candidate mapping information
        labelList mC(cellRemovalList.size(), -1);

        forAll(mC, indexI)
        {
            mC[indexI] = cellRemovalList[indexI];
        }

        // Set the mapping for this cell
        setCellMapping(newCellIndex[cellI], mC);
    }

    // Set fill-in mapping for two new boundary faces
    if (edgePatch > -1)
    {
        forAll(newBdyFaceIndex, i)
        {
            // Set the mapping for this face
            setFaceMapping(newBdyFaceIndex[i], oldBdyFaceIndex);
        }
    }

    // Fill in mapping information for the new face.
    // Since it is internal, interpolate fluxes by default.
    setFaceMapping(newFaceIndex);

    if (debug > 2)
    {
        if (edgePatch > -1)
        {
            Info << "Added edge: " << endl;

            Info << newEdgeIndex << ":: "
                 << edges_[newEdgeIndex]
                 << " edgeFaces: "
                 << edgeFaces_[newEdgeIndex]
                 << endl;
        }

        Info << "Added face(s): " << endl;

        Info << newFaceIndex << ":: "
             << faces_[newFaceIndex];

        Info << " faceEdges: "
             << faceEdges_[newFaceIndex]
             << endl;

        if (edgePatch > -1)
        {
            forAll(newBdyFaceIndex, faceI)
            {
                Info << newBdyFaceIndex[faceI] << ":: "
                     << faces_[newBdyFaceIndex[faceI]]
                     << " faceEdges: "
                     << faceEdges_[newBdyFaceIndex[faceI]]
                     << endl;
            }
        }

        Info << "Added cells: " << endl;

        forAll(newCellIndex, cellI)
        {
            Info << newCellIndex[cellI] << ":: "
                 << cells_[newCellIndex[cellI]]
                 << endl;
        }

        if (debug > 3)
        {
            writeVTK
            (
                Foam::name(eIndex)
              + '(' + Foam::name(edgeToCheck[0])
              + ',' + Foam::name(edgeToCheck[1]) + ')'
              + "_afterSwap_"
              + Foam::name(numTriangulations) + "_"
              + Foam::name(triangulationIndex),
                newCellIndex
            );
        }
    }

    // Specify that the swap was successful
    map.type() = 1;

    // Return the changeMap
    return map;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
