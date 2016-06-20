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

#include "dynamicTopoFvMesh.H"

#include "triFace.H"
#include "objectMap.H"
#include "changeMap.H"
#include "triPointRef.H"
#include "linePointRef.H"
#include "multiThreader.H"
#include "coupledInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Perform a Delaunay test on an internal face
// - Returns 'true' if the test failed
bool dynamicTopoFvMesh::testDelaunay
(
    const label fIndex
) const
{
    bool failed = false, procCoupled = false;
    label eIndex = -1, pointIndex = -1, fLabel = -1;
    label sIndex = -1, pIndex = -1;
    FixedList<triFace, 2> triFaces(triFace(-1, -1, -1));
    FixedList<bool, 2> foundTriFace(false);

    // If this entity was deleted, skip it.
    if (faces_[fIndex].empty())
    {
        return failed;
    }

    // If face is not shared by prism cells, skip it.
    label own = owner_[fIndex], nei = neighbour_[fIndex];

    if (cells_[own].size() != 5)
    {
        return failed;
    }

    if (nei > -1)
    {
        if (cells_[nei].size() != 5)
        {
            return failed;
        }
    }

    // Boundary faces are discarded.
    if (whichPatch(fIndex) > -1)
    {
        procCoupled = processorCoupledEntity(fIndex);

        if (!procCoupled)
        {
            return failed;
        }
    }

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
                }
                else
                {
                    // Update the first face.
                    triFaces[0][0] = thisFace[0];
                    triFaces[0][1] = thisFace[1];
                    triFaces[0][2] = thisFace[2];

                    foundTriFace[0] = true;
                    fLabel = eFaces[faceI];

                    // Take this edge
                    eIndex = fEdges[edgeI];

                    if (procCoupled)
                    {
                        // Stop searching for processor boundary cases
                        foundTriFace[1] = true;
                    }
                }
            }
        }
    }

    const edge& checkEdge = edges_[eIndex];

    // Configure the comparison edge
    edge cEdge(-1, -1);

    if (procCoupled)
    {
        // If this is a new entity, bail out for now.
        // This will be handled at the next time-step.
        if (fIndex >= nOldFaces_)
        {
            return failed;
        }

        const label faceEnum = coupleMap::FACE;
        const label pointEnum = coupleMap::POINT;

        // Check slaves
        bool foundEdge = false;

        forAll(procIndices_, pI)
        {
            // Fetch non-const reference to subMeshes
            const coupledMesh& recvMesh = recvMeshes_[pI];

            const coupleMap& cMap = recvMesh.map();
            const dynamicTopoFvMesh& sMesh = recvMesh.subMesh();

            if ((sIndex = cMap.findSlave(faceEnum, fIndex)) > -1)
            {
                // Find equivalent points on the slave
                cEdge[0] = cMap.findSlave(pointEnum, checkEdge[0]);
                cEdge[1] = cMap.findSlave(pointEnum, checkEdge[1]);

                // Find a triangular face containing cEdge
                const labelList& sfE = sMesh.faceEdges_[sIndex];

                forAll(sfE, edgeI)
                {
                    // Obtain edgeFaces for this edge
                    const labelList& seF = sMesh.edgeFaces_[sfE[edgeI]];

                    forAll(seF, faceI)
                    {
                        const face& tF = sMesh.faces_[seF[faceI]];

                        if (tF.size() == 3)
                        {
                            if
                            (
                                (findIndex(tF, cEdge[0]) > -1) &&
                                (findIndex(tF, cEdge[1]) > -1)
                            )
                            {
                                triFaces[1][0] = tF[0];
                                triFaces[1][1] = tF[1];
                                triFaces[1][2] = tF[2];

                                foundEdge = true;

                                break;
                            }
                        }
                    }

                    if (foundEdge)
                    {
                        break;
                    }
                }

                // Store patch index for later
                pIndex = pI;

                break;
            }
        }

        if (sIndex == -1 || !foundEdge)
        {
            // Write out for post-processing
            writeVTK("coupledDelaunayFace_" + Foam::name(fIndex), fIndex, 2);

            FatalErrorIn
            (
                "bool dynamicTopoFvMesh::testDelaunay"
                "(const label fIndex) const"
            )
                << "Coupled maps were improperly specified." << nl
                << " Slave index not found for: " << nl
                << " Face: " << fIndex << nl
                << " checkEdge: " << checkEdge << nl
                << " cEdge: " << cEdge << nl
                << abort(FatalError);
        }
    }
    else
    {
        // Non-coupled case
        cEdge = checkEdge;
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

    scalar rSquared = (a - cCentre) & (a - cCentre);

    // Find the isolated point on the second face
    point otherPoint = vector::zero;

    // Check the first point
    if (triFaces[1][0] != cEdge[0] && triFaces[1][0] != cEdge[1])
    {
        pointIndex = triFaces[1][0];
    }

    // Check the second point
    if (triFaces[1][1] != cEdge[0] && triFaces[1][1] != cEdge[1])
    {
        pointIndex = triFaces[1][1];
    }

    // Check the third point
    if (triFaces[1][2] != cEdge[0] && triFaces[1][2] != cEdge[1])
    {
        pointIndex = triFaces[1][2];
    }

    if (procCoupled)
    {
        otherPoint = recvMeshes_[pIndex].subMesh().points_[pointIndex];
    }
    else
    {
        otherPoint = points_[pointIndex];
    }

    // ...and determine whether it lies in this circle
    if (((otherPoint - cCentre) & (otherPoint - cCentre)) < rSquared)
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
    FixedList<label,2> c0BdyIndex(-1), c0IntIndex(-1);
    FixedList<label,2> c1BdyIndex(-1), c1IntIndex(-1);
    FixedList<face,2> c0BdyFace(face(3)), c0IntFace(face(4));
    FixedList<face,2> c1BdyFace(face(3)), c1IntFace(face(4));
    FixedList<label,2> commonEdgeIndex(-1);
    FixedList<edge,2> commonEdges(edge(-1, -1));
    FixedList<label,4> otherEdgeIndex(-1);
    FixedList<label,4> commonFaceIndex(-1), cornerEdgeIndex(-1);
    FixedList<face,4> commonFaces(face(3)), commonIntFaces(face(4));
    FixedList<label,4> commonIntFaceIndex(-1);
    FixedList<bool,2> foundTriFace0(false), foundTriFace1(false);
    FixedList<face,2> triFaces0(face(3)), triFaces1(face(3));

    if (coupledModification_)
    {
        if (processorCoupledEntity(fIndex))
        {
            // When swapping across processors, we need to add the
            // prism-cell from (as well as delete on) the patchSubMesh
            const changeMap faceMap = insertCells(fIndex);

            // Bail out if entity is handled elsewhere
            if (faceMap.type() == -2)
            {
                return faceMap;
            }

            if (faceMap.type() != 1)
            {
                FatalErrorIn
                (
                    "const changeMap dynamicTopoFvMesh::swapQuadFace"
                    "(const label fIndex)"
                )
                    << " Could not insert cells around face: " << fIndex
                    << " :: " << faces_[fIndex] << nl
                    << abort(FatalError);
            }

            // Figure out the new internal face index.
            // This should not be a coupled face anymore.
            label newIndex = faceMap.index();

            // Temporarily turn off coupledModification
            unsetCoupledModification();

            // Recursively call this function for the new face
            map = swapQuadFace(newIndex);

            // Turn it back on.
            setCoupledModification();

            return map;
        }
    }

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
        Pout<< nl << nl << "Face: " << fIndex
            << " needs to be flipped. " << endl;

        if (debug > 2)
        {
            Pout<< " Cell[0]: " << c0 << ": " << cells_[c0] << nl
                << " Cell[1]: " << c1 << ": " << cells_[c1] << nl;

            Pout<< " Common Faces: Set 1: "
                << commonFaceIndex[0] << ": " << commonFaces[0] << ", "
                << commonFaceIndex[1] << ": " << commonFaces[1] << nl;

            Pout<< " Common Faces: Set 2: "
                << commonFaceIndex[2] << ": " << commonFaces[2] << ", "
                << commonFaceIndex[3] << ": " << commonFaces[3] << nl;

            Pout<< " Old face: " << faces_[fIndex] << nl
                << " Old faceEdges: " << faceEdges_[fIndex] << endl;
        }

        // Write out VTK files before change
        if (debug > 3)
        {
            labelList cellHull(2, -1);

            cellHull[0] = c0;
            cellHull[1] = c1;

            writeVTK(Foam::name(fIndex) + "_Swap_0", cellHull);
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
        Pout<< "New flipped face: " << newFace << nl;

        if (debug > 2)
        {
            forAll(newBdyFace, faceI)
            {
                Pout<< "New boundary face[" << faceI << "]: "
                    << commonFaceIndex[faceI]
                    << ": " << newBdyFace[faceI] << nl;
            }
        }

        Pout<< endl;
    }

    // Check the orientation of the two quad faces, and modify as necessary
    label newOwn=0, newNei=0;

    // The quad face belonging to cell[1] now becomes a part of cell[0]
    if (neighbour_[commonIntFaceIndex[1]] == -1)
    {
        // Boundary face
        // Face doesn't need to be flipped, just update the owner
        f = commonIntFaces[1];
        newOwn = c0;
        newNei = -1;
    }
    else
    if (owner_[commonIntFaceIndex[1]] == c1)
    {
        // This face is on the interior, check for previous owner
        // Upper-triangular ordering has to be maintained, however...
        if (c0 > neighbour_[commonIntFaceIndex[1]])
        {
            // Flip is necessary
            f = commonIntFaces[1].reverseFace();
            newOwn = neighbour_[commonIntFaceIndex[1]];
            newNei = c0;

            setFlip(commonIntFaceIndex[1]);
        }
        else
        {
            // Flip isn't necessary, just change the owner
            f = commonIntFaces[1];
            newOwn = c0;
            newNei = neighbour_[commonIntFaceIndex[1]];
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
            f = commonIntFaces[1].reverseFace();
            newOwn = c0;
            newNei = owner_[commonIntFaceIndex[1]];

            setFlip(commonIntFaceIndex[1]);
        }
        else
        {
            // Flip isn't necessary, just change the neighbour
            f = commonIntFaces[1];
            newOwn = owner_[commonIntFaceIndex[1]];
            newNei = c0;
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
        f = commonIntFaces[2];
        newOwn = c1;
        newNei = -1;
    }
    else
    if (owner_[commonIntFaceIndex[2]] == c0)
    {
        // This face is on the interior, check for previous owner
        // Upper-triangular ordering has to be maintained, however...
        if (c1 > neighbour_[commonIntFaceIndex[2]])
        {
            // Flip is necessary
            f = commonIntFaces[2].reverseFace();
            newOwn = neighbour_[commonIntFaceIndex[2]];
            newNei = c1;

            setFlip(commonIntFaceIndex[2]);
        }
        else
        {
            // Flip isn't necessary, just change the owner
            f = commonIntFaces[2];
            newOwn = c1;
            newNei = neighbour_[commonIntFaceIndex[2]];
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
            f = commonIntFaces[2].reverseFace();
            newOwn = c1;
            newNei = owner_[commonIntFaceIndex[2]];

            setFlip(commonIntFaceIndex[2]);
        }
        else
        {
            // Flip isn't necessary, just change the neighbour
            f = commonIntFaces[2];
            newOwn = owner_[commonIntFaceIndex[2]];
            newNei = c1;
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
    status(TOTAL_SWAPS)++;

    // Return a successful operation.
    map.type() = 1;

    return map;
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
    label numIndices = coupledModification_ ? 2 : 1;

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
    scalar& minQuality,
    labelList& m,
    labelList& hullVertices,
    PtrList<scalarListList>& Q,
    PtrList<labelListList>& K,
    PtrList<labelListList>& triangulations,
    const label checkIndex
) const
{
    // Obtain a reference to this edge
    const labelList& edgeFaces = edgeFaces_[eIndex];

    // If this entity was deleted, skip it.
    if (edgeFaces.empty())
    {
        return false;
    }

    // Ensure that edge is surrounded by triangles
    forAll(edgeFaces, faceI)
    {
        if (faces_[edgeFaces[faceI]].size() != 3)
        {
            return false;
        }
    }

    const edge& edgeToCheck = edges_[eIndex];

    if (coupledModification_)
    {
        return coupledFillTables(eIndex, minQuality, m, Q, K, triangulations);
    }

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

    // Fill dynamic programming tables
    fillTables
    (
        edgeToCheck,
        minQuality,
        m[checkIndex],
        hullVertices,
        points_,
        Q[checkIndex],
        K[checkIndex],
        triangulations[checkIndex]
    );

    // Print out tables for debugging
    if (debug > 5)
    {
        printTables(m, Q, K, checkIndex);
    }

    return true;
}


// Fill tables given addressing
void dynamicTopoFvMesh::fillTables
(
    const edge& edgeToCheck,
    const scalar minQuality,
    const label m,
    const labelList& hullVertices,
    const UList<point>& points,
    scalarListList& Q,
    labelListList& K,
    labelListList& triangulations
) const
{
    for (label i = (m - 3); i >= 0; i--)
    {
        for (label j = i + 2; j < m; j++)
        {
            for (label k = i + 1; k < j; k++)
            {
                scalar q = (*tetMetric_)
                (
                    points[hullVertices[i]],
                    points[hullVertices[k]],
                    points[hullVertices[j]],
                    points[edgeToCheck[0]]
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
                                points[hullVertices[j]],
                                points[hullVertices[k]],
                                points[hullVertices[i]],
                                points[edgeToCheck[1]]
                            )
                        )
                    );
                }

                if (k < j - 1)
                {
                    q = Foam::min(q, Q[k][j]);
                }

                if (k > i + 1)
                {
                    q = Foam::min(q, Q[i][k]);
                }

                if ((k == i + 1) || (q > Q[i][j]))
                {
                    Q[i][j] = q;
                    K[i][j] = k;
                }
            }
        }
    }
}


// Remove the edge according to the swap sequence.
// - Returns a changeMap with a type specifying:
//     1: Swap sequence was successful
//    -1: Swap sequence failed
const changeMap dynamicTopoFvMesh::removeEdgeFlips
(
    const label eIndex,
    const scalar minQuality,
    const labelList& vertexHull,
    PtrList<scalarListList>& Q,
    PtrList<labelListList>& K,
    PtrList<labelListList>& triangulations,
    const label checkIndex
)
{
    changeMap map, slaveMap;

    if (debug > 2)
    {
        Pout<< nl << " Removing edge : " << eIndex << " by flipping."
            << " Edge: " << edges_[eIndex]
            << " minQuality: " << minQuality
            << endl;
    }

    if (coupledModification_)
    {
        if (processorCoupledEntity(eIndex))
        {
            const changeMap edgeMap = insertCells(eIndex);

            // Bail out if entity is handled elsewhere
            if (edgeMap.type() == -2)
            {
                return edgeMap;
            }

            if (edgeMap.type() != 1)
            {
                FatalErrorIn
                (
                    "\n"
                    "const changeMap dynamicTopoFvMesh::removeEdgeFlips\n"
                    "(\n"
                    "    const label eIndex,\n"
                    "    const scalar minQuality,\n"
                    "    const labelList& vertexHull,\n"
                    "    PtrList<scalarListList>& Q,\n"
                    "    PtrList<labelListList>& K,\n"
                    "    PtrList<labelListList>& triangulations,\n"
                    "    const label checkIndex\n"
                    ")\n"
                )
                    << " Could not insert cells around edge: " << eIndex
                    << " :: " << edges_[eIndex] << nl
                    << abort(FatalError);
            }

            // Temporarily turn off coupledModification
            unsetCoupledModification();

            // Re-fill tables for the reconfigured edge.
            // This should not be a coupled edge anymore.
            label newIndex = edgeMap.index();

            if (processorCoupledEntity(newIndex))
            {
                // Write out edge connectivity
                writeEdgeConnectivity(newIndex);

                FatalErrorIn
                (
                    "\n"
                    "const changeMap dynamicTopoFvMesh::removeEdgeFlips\n"
                    "(\n"
                    "    const label eIndex,\n"
                    "    const scalar minQuality,\n"
                    "    const labelList& vertexHull,\n"
                    "    PtrList<scalarListList>& Q,\n"
                    "    PtrList<labelListList>& K,\n"
                    "    PtrList<labelListList>& triangulations,\n"
                    "    const label checkIndex\n"
                    ")\n"
                )
                    << " Edge: " << newIndex
                    << " :: " << edges_[newIndex] << nl
                    << " is still processor-coupled. "
                    << abort(FatalError);
            }

            const edge& newEdge = edges_[newIndex];

            // Build vertexHull for this edge
            labelList newVertexHull;
            buildVertexHull(newIndex, newVertexHull);

            fillTables
            (
                newEdge,
                minQuality,
                newVertexHull.size(),
                newVertexHull,
                points_,
                Q[0],
                K[0],
                triangulations[0]
            );

            // Recursively call this function for the new edge
            map =
            (
                removeEdgeFlips
                (
                    newIndex,
                    minQuality,
                    newVertexHull,
                    Q,
                    K,
                    triangulations
                )
            );

            // Turn it back on.
            setCoupledModification();

            return map;
        }
    }

    label m = vertexHull.size();
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
        vertexHull,
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
            vertexHull,
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
            vertexHull,
            triangulations[checkIndex]
        )
    );

    // Check that the triangulation is valid
    label pIndex = -1;

    if (tF == -1)
    {
        const edge& edgeToCheck = edges_[eIndex];

        Pout<< " All triangulations: " << nl
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
            "    const labelList& vertexHull,\n"
            "    PtrList<scalarListList>& Q,\n"
            "    PtrList<labelListList>& K,\n"
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
        Pout<< " Identified tF as: " << tF << nl;

        Pout<< " Triangulation: "
            << triangulations[checkIndex][0][tF] << " "
            << triangulations[checkIndex][1][tF] << " "
            << triangulations[checkIndex][2][tF] << " "
            << nl;

        Pout<< " All triangulations: " << nl
            << ' ' << triangulations[checkIndex][0] << nl
            << ' ' << triangulations[checkIndex][1] << nl
            << ' ' << triangulations[checkIndex][2] << nl
            << endl;
    }

    if (coupledModification_)
    {
        if (locallyCoupledEntity(eIndex))
        {
            // Flip the slave edge as well.
            label sIndex = -1;

            // Determine the slave index.
            forAll(patchCoupling_, patchI)
            {
                if (patchCoupling_(patchI))
                {
                    const label edgeEnum  = coupleMap::EDGE;
                    const coupleMap& cMap = patchCoupling_[patchI].map();

                    if ((sIndex = cMap.findSlave(edgeEnum, eIndex)) > -1)
                    {
                        break;
                    }
                }
            }

            if (sIndex == -1)
            {
                FatalErrorIn
                (
                    "\n"
                    "const changeMap dynamicTopoFvMesh::removeEdgeFlips\n"
                    "(\n"
                    "    const label eIndex,\n"
                    "    const scalar minQuality,\n"
                    "    const labelList& vertexHull,\n"
                    "    PtrList<scalarListList>& Q,\n"
                    "    PtrList<labelListList>& K,\n"
                    "    PtrList<labelListList>& triangulations,\n"
                    "    const label checkIndex\n"
                    ")\n"
                )
                    << "Coupled maps were improperly specified." << nl
                    << " Slave index not found for: " << nl
                    << " Edge: " << eIndex << nl
                    << abort(FatalError);
            }

            if (debug > 2)
            {
                Pout<< nl << "Removing slave edge: " << sIndex
                    << " for master edge: " << eIndex << endl;
            }

            // Build vertexHull for this edge
            labelList slaveVertexHull;
            buildVertexHull(sIndex, slaveVertexHull);

            // Turn off switch temporarily.
            unsetCoupledModification();

            // Recursively call for the slave edge.
            slaveMap =
            (
                removeEdgeFlips
                (
                    sIndex,
                    minQuality,
                    slaveVertexHull,
                    Q,
                    K,
                    triangulations,
                    1
                )
            );

            // Turn it back on.
            setCoupledModification();

            // Bail out if the slave failed.
            if (slaveMap.type() == -1)
            {
                // Reset all triangulations and bail out
                triangulations[checkIndex][0] = -1;
                triangulations[checkIndex][1] = -1;
                triangulations[checkIndex][2] = -1;

                return slaveMap;
            }
        }
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
                            vertexHull,
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
            Pout<< "Triangulations: " << nl;

            forAll(triangulations[checkIndex], row)
            {
                Pout<< triangulations[checkIndex][row] << nl;
            }

            Pout<<endl;

            // Should have performed at least one swap
            FatalErrorIn
            (
                "\n"
                "const changeMap dynamicTopoFvMesh::removeEdgeFlips\n"
                "(\n"
                "    const label eIndex,\n"
                "    const scalar minQuality,\n"
                "    const labelList& vertexHull,\n"
                "    PtrList<scalarListList>& Q,\n"
                "    PtrList<labelListList>& K,\n"
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
            vertexHull,
            hullFaces,
            hullCells
        )
    );

    // Done with this face, so reset it
    triangulations[checkIndex][0][tF] = -1;
    triangulations[checkIndex][1][tF] = -1;
    triangulations[checkIndex][2][tF] = -1;

    // Update the coupled map
    if (coupledModification_)
    {
        // Create a mapping entry for the new edge.
        const coupleMap& cMap = patchCoupling_[pIndex].map();

        if (locallyCoupledEntity(map.addedEdgeList()[0].index()))
        {
            cMap.mapSlave
            (
                coupleMap::EDGE,
                map.addedEdgeList()[0].index(),
                slaveMap.addedEdgeList()[0].index()
            );

            cMap.mapMaster
            (
                coupleMap::EDGE,
                slaveMap.addedEdgeList()[0].index(),
                map.addedEdgeList()[0].index()
            );
        }

        // Add a mapping entry for two new faces as well.
        face cF(3);

        const List<objectMap>& amfList = map.addedFaceList();
        const List<objectMap>& asfList = slaveMap.addedFaceList();

        forAll(amfList, mfI)
        {
            // Configure a face for comparison.
            const face& mF = faces_[amfList[mfI].index()];

            forAll(mF, pointI)
            {
                cF[pointI] = cMap.entityMap(coupleMap::POINT)[mF[pointI]];
            }

            bool matched = false;

            forAll(asfList, sfI)
            {
                const face& sF = faces_[asfList[sfI].index()];

                if (triFace::compare(triFace(cF), triFace(sF)))
                {
                    cMap.mapSlave
                    (
                        coupleMap::FACE,
                        amfList[mfI].index(),
                        asfList[sfI].index()
                    );

                    cMap.mapMaster
                    (
                        coupleMap::FACE,
                        asfList[sfI].index(),
                        amfList[mfI].index()
                    );

                    matched = true;

                    break;
                }
            }

            if (!matched)
            {
                Pout<< "masterFaces: " << nl
                    << amfList << endl;

                Pout<< "slaveFaces: " << nl
                    << asfList << endl;

                forAll(amfList, mfI)
                {
                    Pout<< amfList[mfI].index() << ": "
                        << faces_[amfList[mfI].index()]
                        << nl;
                }

                forAll(asfList, sfI)
                {
                    Pout<< asfList[sfI].index() << ": "
                        << faces_[asfList[sfI].index()]
                        << nl;
                }

                Pout<< endl;

                FatalErrorIn
                (
                    "\n"
                    "const changeMap dynamicTopoFvMesh::removeEdgeFlips\n"
                    "(\n"
                    "    const label eIndex,\n"
                    "    const scalar minQuality,\n"
                    "    const labelList& vertexHull,\n"
                    "    PtrList<scalarListList>& Q,\n"
                    "    PtrList<labelListList>& K,\n"
                    "    PtrList<labelListList>& triangulations,\n"
                    "    const label checkIndex\n"
                    ")\n"
                )
                    << "Failed to build coupled face maps."
                    << abort(FatalError);
            }
        }
    }

    // Finally remove the edge
    removeEdge(eIndex);

    // Update map
    map.removeEdge(eIndex);

    // Increment the counter
    status(TOTAL_SWAPS)++;

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
    const labelListList& triangulations,
    bool output
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
                intPt
            )
        );

        if (intersects)
        {
            return i;
        }
    }

    if (debug > 1 || output)
    {
        Pout<< nl << nl << "Hull Vertices: " << nl;

        forAll(hullVertices, vertexI)
        {
            Pout<< hullVertices[vertexI] << ": "
                << points_[hullVertices[vertexI]]
                << nl;
        }

        InfoIn
        (
            "\n"
            "label dynamicTopoFvMesh::identify32Swap\n"
            "(\n"
            "    const label eIndex,\n"
            "    const labelList& hullVertices,\n"
            "    const labelListList& triangulations,\n"
            "    bool output\n"
            ") const\n"
        )   << nl
            << "Could not determine 3-2 swap triangulation." << nl
            << "Edge: " << edgeToCheck << nl
            << "Edge Points: "
            << segment.start() << "," << segment.end() << nl
            << endl;
    }

    // Could not find an intersecting triangulation.
    //  - If this is a boundary edge, a curved surface-mesh
    //    was probably the reason why this failed.
    //  - If so, declare the nearest triangulation instead.
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

    if (debug > 1 || output)
    {
        Pout<< " All distances :" << dist << nl
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
    const label eIndex,
    labelList& hullVertices
) const
{
    scalar minQuality = GREAT;

    // Obtain a reference to this edge
    const edge& edgeToCheck = edges_[eIndex];
    const labelList& edgeFaces = edgeFaces_[eIndex];

    // If this entity was deleted, skip it.
    if (edgeFaces.empty())
    {
        return minQuality;
    }

    // Ensure that edge is surrounded by triangles
    forAll(edgeFaces, faceI)
    {
        if (faces_[edgeFaces[faceI]].size() != 3)
        {
            return minQuality;
        }
    }

    // Build vertexHull for this edge
    buildVertexHull(eIndex, hullVertices);

    if (coupledModification_)
    {
        if (locallyCoupledEntity(eIndex))
        {
            // Compute the minimum quality of the slave edge as well.
            label sIndex = -1;

            // Determine the slave index.
            forAll(patchCoupling_, patchI)
            {
                if (patchCoupling_(patchI))
                {
                    const label edgeEnum  = coupleMap::EDGE;
                    const coupleMap& cMap = patchCoupling_[patchI].map();

                    if ((sIndex = cMap.findSlave(edgeEnum, eIndex)) > -1)
                    {
                        break;
                    }
                }
            }

            if (sIndex == -1)
            {
                FatalErrorIn
                (
                    "scalar dynamicTopoFvMesh::computeMinQuality"
                    "(const label eIndex, labelList& hullVertices) const"
                )
                    << nl << "Coupled maps were improperly specified." << nl
                    << " Slave index not found for: " << nl
                    << " Edge: " << eIndex << nl
                    << abort(FatalError);
            }

            // Build vertexHull for this edge
            labelList slaveVertexHull;
            buildVertexHull(eIndex, slaveVertexHull);

            // Temporarily turn off coupledModification
            unsetCoupledModification();

            scalar slaveQuality = computeMinQuality(sIndex, slaveVertexHull);

            minQuality = Foam::min(slaveQuality, minQuality);

            // Turn it back on.
            setCoupledModification();
        }
        else
        if (processorCoupledEntity(eIndex))
        {
            // Don't compute minQuality here, but do it in
            // fillTables instead, while calculating hullPoints
            return minQuality;
        }
    }

    // Compute minQuality
    minQuality =
    (
        Foam::min
        (
            computeMinQuality
            (
                edgeToCheck,
                hullVertices,
                points_,
                (whichEdgePatch(eIndex) < 0)
            ),
            minQuality
        )
    );

    // Ensure that the mesh is valid
    if (minQuality < 0.0)
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

        FatalErrorIn
        (
            "scalar dynamicTopoFvMesh::computeMinQuality"
            "(const label eIndex, labelList& hullVertices) const"
        )
            << "Encountered negative cell-quality!" << nl
            << "Edge: " << eIndex << ": " << edgeToCheck << nl
            << "vertexHull: " << hullVertices << nl
            << "Minimum Quality: " << minQuality
            << abort(FatalError);
    }

    return minQuality;
}


// Compute minQuality given addressing
scalar dynamicTopoFvMesh::computeMinQuality
(
    const edge& edgeToCheck,
    const labelList& hullVertices,
    const UList<point>& points,
    bool closedRing
) const
{
    scalar cQuality = 0.0;
    scalar minQuality = GREAT;

    // Obtain point references
    const point& a = points[edgeToCheck[0]];
    const point& c = points[edgeToCheck[1]];

    label start = (closedRing ? 0 : 1);

    for (label indexJ = start; indexJ < hullVertices.size(); indexJ++)
    {
        label indexI = hullVertices.rcIndex(indexJ);

        // Pick vertices off the list
        const point& b = points[hullVertices[indexI]];
        const point& d = points[hullVertices[indexJ]];

        // Compute the quality
        cQuality = tetMetric_(a, b, c, d);

        // Check if the quality is worse
        minQuality = Foam::min(cQuality, minQuality);
    }

    return minQuality;
}


// Method used to perform a 2-3 swap in 3D
// - Returns a changeMap with a type specifying:
//     1: Swap was successful
// - The index of the triangulated face in map.index()
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
    //      [1] Remove face: [ edgeToCheck[0] edgeToCheck[1] isolatedVertex ]
    //      [2] Remove two cells on either side of removed face
    //      [3] Add one edge
    //      [4] Add three new faces
    //      [5] Add three new cells
    //      Update faceEdges and edgeFaces information

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
        Pout<< nl
            << "== Swapping 2-3 ==" << nl
            << "Edge: " << eIndex << ": " << edgeToCheck << endl;

        if (debug > 2)
        {
            Pout<< " On SubMesh: " << isSubMesh_ << nl;
            Pout<< " coupledModification: " << coupledModification_ << nl;

            label bPatch = whichEdgePatch(eIndex);

            const polyBoundaryMesh& boundary = boundaryMesh();

            if (bPatch == -1)
            {
                Pout<< " Patch: Internal" << nl;
            }
            else
            if (bPatch < boundary.size())
            {
                Pout<< " Patch: " << boundary[bPatch].name() << nl;
            }
            else
            {
                Pout<< " New patch: " << bPatch << endl;
            }

            Pout<< " Ring: " << hullVertices << nl
                << " Faces: " << hullFaces << nl
                << " Cells: " << hullCells << nl
                << " Triangulation: "
                << triangulations[0][triangulationIndex] << " "
                << triangulations[1][triangulationIndex] << " "
                << triangulations[2][triangulationIndex] << " "
                << nl
                << " Isolated vertex: " << isolatedVertex << endl;
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
                labelList(cellsForRemoval),
                3, false, true
            );
        }
    }

    // Check if this is an internal face
    if (cellsForRemoval[1] == -1)
    {
        // Write out for post-processing
        Pout<< " isolatedVertex: " << isolatedVertex << nl
            << " triangulations: " << triangulations << nl
            << " numTriangulations: " << numTriangulations << nl
            << " triangulationIndex: " << triangulationIndex << endl;

        writeVTK("Edge23_" + Foam::name(eIndex), eIndex, 1);
        writeVTK("Cells23_" + Foam::name(eIndex), hullCells, 3);

        // Write out identify32Swap output for diagnostics
        identify32Swap(eIndex, hullVertices, triangulations, true);

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

        // Add this cell to the map.
        map.addCell(newCellIndex[cellI]);
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
            newCellIndex[1],
            labelList(0)
        )
    );

    // Add this face to the map.
    map.addFace(newFaceIndex[0]);

    // Note the triangulation face in index()
    map.index() = newFaceIndex[0];

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
            newCellIndex[2],
            labelList(0)
        )
    );

    // Add this face to the map.
    map.addFace(newFaceIndex[1]);

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
            newCellIndex[2],
            labelList(0)
        )
    );

    // Add this face to the map.
    map.addFace(newFaceIndex[2]);

    // Add an entry to edgeFaces
    labelList newEdgeFaces(3);
    newEdgeFaces[0] = newFaceIndex[0];
    newEdgeFaces[1] = newFaceIndex[1];
    newEdgeFaces[2] = newFaceIndex[2];

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
            newEdgeFaces
        )
    );

    // Add this edge to the map.
    map.addEdge(newEdgeIndex);

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

                // Update faceEdges and edgeFaces
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
                    }

                    if (edges_[faceEdges_[faceIndex][edgeI]] == check[1])
                    {
                        newFaceEdges[0][nE0++] = faceEdges_[faceIndex][edgeI];

                        meshOps::sizeUpList
                        (
                            newFaceIndex[0],
                            edgeFaces_[faceEdges_[faceIndex][edgeI]]
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
                    }

                    if (edges_[faceEdges_[faceIndex][edgeI]] == check[4])
                    {
                        newFaceEdges[1][nE1++] = faceEdges_[faceIndex][edgeI];

                        meshOps::sizeUpList
                        (
                            newFaceIndex[1],
                            edgeFaces_[faceEdges_[faceIndex][edgeI]]
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

                // Update faceEdges and edgeFaces
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
                    }

                    if (edges_[fEdges[edgeI]] == check[5])
                    {
                        newFaceEdges[2][nE2++] = fEdges[edgeI];

                        meshOps::sizeUpList
                        (
                            newFaceIndex[2],
                            edgeFaces_[fEdges[edgeI]]
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

    // Update edgeFaces for edges of the removed face
    forAll(faceEdges_[faceForRemoval], edgeI)
    {
        label edgeIndex = faceEdges_[faceForRemoval][edgeI];

        meshOps::sizeDownList
        (
            faceForRemoval,
            edgeFaces_[edgeIndex]
        );
    }

    // Remove the face
    removeFace(faceForRemoval);

    // Update map
    map.removeFace(faceForRemoval);

    forAll(cellsForRemoval, cellI)
    {
        removeCell(cellsForRemoval[cellI]);

        // Update map
        map.removeCell(cellsForRemoval[cellI]);
    }

    // Update the cell list with newly configured cells.
    forAll(newCellIndex, cellI)
    {
        cells_[newCellIndex[cellI]] = newTetCell[cellI];

        if (cellI == 2)
        {
            // Skip mapping for the intermediate cell.
            setCellMapping(newCellIndex[cellI], hullCells, false);
        }
        else
        {
            // Set the mapping for this cell
            setCellMapping(newCellIndex[cellI], hullCells);
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
        Pout<< "Added edge: " << nl;

        Pout<< newEdgeIndex << ":: "
            << edges_[newEdgeIndex]
            << " edgeFaces: "
            << edgeFaces_[newEdgeIndex]
            << nl;

        Pout<< "Added faces: " << nl;

        forAll(newFaceIndex, faceI)
        {
            Pout<< newFaceIndex[faceI] << ":: "
                << faces_[newFaceIndex[faceI]]
                << " faceEdges: "
                << faceEdges_[newFaceIndex[faceI]]
                << nl;
        }

        Pout<< "Added cells: " << nl;

        forAll(newCellIndex, cellI)
        {
            Pout<< newCellIndex[cellI] << ":: "
                << cells_[newCellIndex[cellI]]
                << nl;
        }

        Pout<< endl;

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
                labelList(newCellIndex),
                3, false, true
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
// - The index of the triangulated face in map.index()
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
    //      Update faceEdges and edgeFaces information

    changeMap map;

    // Obtain a copy of the edge
    edge edgeToCheck = edges_[eIndex];

    // Determine the patch this edge belongs to
    label edgePatch = whichEdgePatch(eIndex);

    // Determine the three faces to be removed
    FixedList<label,3> facesForRemoval;
    dynamicLabelList cellRemovalList(3);

    forAll(facesForRemoval, faceI)
    {
        facesForRemoval[faceI] =
        (
            hullFaces[triangulations[faceI][triangulationIndex]]
        );

        label own = owner_[facesForRemoval[faceI]];
        label nei = neighbour_[facesForRemoval[faceI]];

        // Check and add cells as well
        if (findIndex(cellRemovalList, own) == -1)
        {
            cellRemovalList.append(own);
        }

        if (nei != -1)
        {
            if (findIndex(cellRemovalList, nei) == -1)
            {
                cellRemovalList.append(nei);
            }
        }
    }

    if (debug > 1)
    {
        // Print out arguments
        Pout<< nl;

        if (edgePatch < 0)
        {
            Pout<< "== Swapping 3-2 ==" << nl;
        }
        else
        {
            Pout<< "== Swapping 2-2 ==" << nl;
        }

        Pout<< " Edge: " << eIndex << ": " << edgeToCheck << endl;

        if (debug > 2)
        {
            Pout<< " On SubMesh: " << isSubMesh_ << nl;
            Pout<< " coupledModification: " << coupledModification_ << nl;

            label bPatch = whichEdgePatch(eIndex);

            const polyBoundaryMesh& boundary = boundaryMesh();

            if (bPatch == -1)
            {
                Pout<< " Patch: Internal" << nl;
            }
            else
            if (bPatch < boundary.size())
            {
                Pout<< " Patch: " << boundary[bPatch].name() << nl;
            }
            else
            {
                Pout<< " New patch: " << bPatch << endl;
            }

            Pout<< " Ring: " << hullVertices << nl
                << " Faces: " << hullFaces << nl
                << " Cells: " << hullCells << nl
                << " Triangulation: "
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
                cellRemovalList,
                3, false, true
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

        // Add this cell to the map.
        map.addCell(newCellIndex[cellI]);
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
            newCellIndex[1],
            labelList(3, -1)
        )
    );

    // Add this face to the map.
    map.addFace(newFaceIndex);

    // Note the triangulation face in index()
    map.index() = newFaceIndex;

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
                -1,
                labelList(3, -1)
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
                -1,
                labelList(3, -1)
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

        // Insert the edge
        newEdgeIndex =
        (
            insertEdge
            (
                edgePatch,
                newEdge,
                newBdyEdgeFaces
            )
        );

        // Add this edge to the map.
        map.addEdge(newEdgeIndex);

        // Update faceEdges with the new edge
        faceEdges_[newFaceIndex][nE++] = newEdgeIndex;
        bdyFaceEdges[0][nBE[0]++] = newEdgeIndex;
        bdyFaceEdges[1][nBE[1]++] = newEdgeIndex;

        // Update edgeFaces with the two new faces
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
        }

        // Add faceEdges for the two new boundary faces
        faceEdges_[newBdyFaceIndex[0]] = bdyFaceEdges[0];
        faceEdges_[newBdyFaceIndex[1]] = bdyFaceEdges[1];

        // Update the number of surface swaps.
        status(SURFACE_SWAPS)++;
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

                // Update faceEdges and edgeFaces
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
        // Update edgeFaces
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
            }
        }

        // Now remove the face...
        removeFace(facesForRemoval[faceI]);

        // Update map
        map.removeFace(facesForRemoval[faceI]);
    }

    forAll(cellRemovalList, cellI)
    {
        removeCell(cellRemovalList[cellI]);

        // Update map
        map.removeCell(cellRemovalList[cellI]);
    }

    // Update the cell list with newly configured cells.
    forAll(newCellIndex, cellI)
    {
        cells_[newCellIndex[cellI]] = newTetCell[cellI];

        // Set the mapping for this cell
        setCellMapping(newCellIndex[cellI], hullCells);
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
            Pout<< "Added edge: "
                << newEdgeIndex << ":: "
                << edges_[newEdgeIndex]
                << " edgeFaces: "
                << edgeFaces_[newEdgeIndex]
                << endl;
        }

        Pout<< "Added face(s): " << nl;

        Pout<< newFaceIndex << ":: "
            << faces_[newFaceIndex];

        Pout<< " faceEdges: "
            << faceEdges_[newFaceIndex]
            << endl;

        if (edgePatch > -1)
        {
            forAll(newBdyFaceIndex, faceI)
            {
                Pout<< newBdyFaceIndex[faceI] << ":: "
                    << faces_[newBdyFaceIndex[faceI]]
                    << " faceEdges: "
                    << faceEdges_[newBdyFaceIndex[faceI]]
                    << nl;
            }
        }

        Pout<< "Added cells: " << nl;

        forAll(newCellIndex, cellI)
        {
            Pout<< newCellIndex[cellI] << ":: "
                << cells_[newCellIndex[cellI]]
                << nl;
        }

        Pout<< endl;

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
                labelList(newCellIndex),
                3, false, true
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
