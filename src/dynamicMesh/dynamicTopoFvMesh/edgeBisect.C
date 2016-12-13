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

#include "Stack.H"
#include "triFace.H"
#include "objectMap.H"
#include "changeMap.H"
#include "coupledInfo.H"
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
//    -3: Bisection failed since edge was on a noRefinement patch.
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
    changeMap map;
    List<changeMap> slaveMaps;
    bool bisectingSlave = false;

    if
    (
        (status(TOTAL_MODIFICATIONS) > maxModifications_) &&
        (maxModifications_ > -1)
    )
    {
        // Reached the max allowable topo-changes.
        stack(tIndex).clear();

        return map;
    }

    // Check if edgeRefinements are to be avoided on patch.
    if (baseMesh_.lengthEstimator().checkRefinementPatch(whichPatch(fIndex)))
    {
        map.type() = -3;

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
            "    bool checkOnly,\n"
            "    bool forceOp\n"
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
    FixedList<edge,4> commonEdges(edge(-1, -1));
    FixedList<label,4> cornerEdgeIndex(-1), commonEdgeIndex(-1);
    FixedList<label,4> commonFaceIndex(-1);
    FixedList<label,2> newPointIndex(-1), newCellIndex(-1);
    FixedList<label,4> otherEdgeIndex(-1), otherEdgePoint(-1);
    FixedList<label,4> otherPointIndex(-1), nextToOtherPoint(-1);
    FixedList<label,2> c0BdyIndex(-1), c0IntIndex(-1);
    FixedList<label,2> c1BdyIndex(-1), c1IntIndex(-1);
    FixedList<face,2> c0BdyFace(face(3)), c0IntFace(face(4));
    FixedList<face,2> c1BdyFace(face(3)), c1IntFace(face(4));

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

    // If coupled modification is set, and this is a
    // master edge, bisect its slaves as well.
    bool localCouple = false, procCouple = false;

    if (coupledModification_)
    {
        const label faceEnum = coupleMap::FACE;
        const label pointEnum = coupleMap::POINT;

        // Is this a locally coupled face (either master or slave)?
        if (locallyCoupledEntity(fIndex, true))
        {
            localCouple = true;
            procCouple = false;
        }
        else
        if (processorCoupledEntity(fIndex))
        {
            procCouple = true;
            localCouple = false;
        }

        if (localCouple && !procCouple)
        {
            // Determine the slave index.
            label sIndex = -1, pIndex = -1;

            forAll(patchCoupling_, patchI)
            {
                if (patchCoupling_(patchI))
                {
                    const coupleMap& cMap = patchCoupling_[patchI].map();

                    if ((sIndex = cMap.findSlave(faceEnum, fIndex)) > -1)
                    {
                        pIndex = patchI;

                        break;
                    }

                    // The following bit happens only during the sliver
                    // exudation process, since slave faces are
                    // usually not added to the coupled face-stack.
                    if ((sIndex = cMap.findMaster(faceEnum, fIndex)) > -1)
                    {
                        pIndex = patchI;

                        // Notice that we are collapsing a slave face first.
                        bisectingSlave = true;

                        break;
                    }
                }
            }

            if (sIndex == -1)
            {
                FatalErrorIn
                (
                    "\n"
                    "const changeMap dynamicTopoFvMesh::bisectQuadFace\n"
                    "(\n"
                    "    const label fIndex,\n"
                    "    const changeMap& masterMap,\n"
                    "    bool checkOnly,\n"
                    "    bool forceOp"
                    ")\n"
                )
                    << "Coupled maps were improperly specified." << nl
                    << " Slave index not found for: " << nl
                    << " Face: " << fIndex << ": " << faces_[fIndex] << nl
                    << abort(FatalError);
            }
            else
            {
                // If we've found the slave, size up the list
                meshOps::sizeUpList
                (
                    changeMap(),
                    slaveMaps
                );

                // Save index and patch for posterity
                slaveMaps[0].index() = sIndex;
                slaveMaps[0].patchIndex() = pIndex;
            }

            if (debug > 1)
            {
                Pout<< nl << " >> Bisecting slave face: " << sIndex
                    << " for master face: " << fIndex << endl;
            }
        }
        else
        if (procCouple && !localCouple)
        {
            // If this is a new entity, bail out for now.
            // This will be handled at the next time-step.
            if (fIndex >= nOldFaces_)
            {
                return map;
            }

            // Check slaves
            forAll(procIndices_, pI)
            {
                // Fetch reference to subMesh
                const coupledMesh& recvMesh = recvMeshes_[pI];
                const coupleMap& cMap = recvMesh.map();

                label sIndex = -1;

                if ((sIndex = cMap.findSlave(faceEnum, fIndex)) > -1)
                {
                    if (debug > 3)
                    {
                        Pout<< "Checking slave face: " << sIndex
                            << " on proc: " << procIndices_[pI]
                            << " for master face: " << fIndex
                            << endl;
                    }

                    // Check if a lower-ranked processor is
                    // handling this edge
                    if
                    (
                        priority
                        (
                            procIndices_[pI],
                            lessOp<label>(),
                            Pstream::myProcNo()
                        )
                    )
                    {
                        if (debug > 3)
                        {
                            Pout<< "Face: " << fIndex
                                << " is handled by proc: "
                                << procIndices_[pI]
                                << ", so bailing out."
                                << endl;
                        }

                        return map;
                    }

                    label curIndex = slaveMaps.size();

                    // Size up the list
                    meshOps::sizeUpList
                    (
                        changeMap(),
                        slaveMaps
                    );

                    // Save index and patch for posterity
                    slaveMaps[curIndex].index() = sIndex;
                    slaveMaps[curIndex].patchIndex() = pI;

                    // Only one slave coupling is possible, so bail out
                    break;
                }
            }
        }
        else
        {
            // Something's wrong with coupling maps
            FatalErrorIn
            (
                "\n"
                "const changeMap dynamicTopoFvMesh::bisectQuadFace\n"
                "(\n"
                "    const label fIndex,\n"
                "    const changeMap& masterMap,\n"
                "    bool checkOnly,\n"
                "    bool forceOp"
                ")\n"
            )
                << "Coupled maps were improperly specified." << nl
                << " localCouple: " << localCouple << nl
                << " procCouple: " << procCouple << nl
                << " Face: " << fIndex << ": " << faces_[fIndex] << nl
                << abort(FatalError);
        }

        // Alias for convenience...
        changeMap& slaveMap = slaveMaps[0];

        label sIndex = slaveMap.index();
        label pI = slaveMap.patchIndex();
        const coupleMap* cMapPtr = NULL;

        // Temporarily turn off coupledModification
        unsetCoupledModification();

        if (localCouple)
        {
            cMapPtr = &(patchCoupling_[pI].map());

            // First check the slave for bisection feasibility.
            slaveMap = bisectQuadFace(sIndex, changeMap(), true, forceOp);
        }
        else
        if (procCouple)
        {
            cMapPtr = &(recvMeshes_[pI].map());

            coupledMesh& recvMesh = recvMeshes_[pI];

            // First check the slave for bisection feasibility.
            slaveMap =
            (
                recvMesh.subMesh().bisectQuadFace
                (
                    sIndex,
                    changeMap(),
                    true,
                    forceOp
                )
            );
        }

        // Turn it back on.
        setCoupledModification();

        if (slaveMap.type() <= 0)
        {
            // Slave couldn't perform bisection.
            map.type() = -2;

            return map;
        }

        // Save index and patch for posterity
        slaveMap.index() = sIndex;
        slaveMap.patchIndex() = pI;

        // Alias for convenience..
        const coupleMap& cMap = *cMapPtr;

        // Temporarily turn off coupledModification
        unsetCoupledModification();

        // First check the master for bisection feasibility.
        changeMap masterMap = bisectQuadFace(fIndex, changeMap(), true);

        // Turn it back on.
        setCoupledModification();

        // Master couldn't perform bisection
        if (masterMap.type() <= 0)
        {
            return masterMap;
        }

        // Fill the masterMap with points that
        // we seek maps for...
        FixedList<labelList, 2> slaveEdges(labelList(2, -1));

        forAll(slaveEdges, edgeI)
        {
            labelList& slaveEdge = slaveEdges[edgeI];

            // Renumber to slave indices
            forAll(slaveEdge, pointI)
            {
                slaveEdge[pointI] =
                (
                    cMap.findSlave
                    (
                        pointEnum,
                        commonEdges[edgeI][pointI]
                    )
                );
            }

            masterMap.addPoint(-1, slaveEdge);
        }

        // Temporarily turn off coupledModification
        unsetCoupledModification();

        if (localCouple)
        {
            // Bisect the local slave.
            slaveMap = bisectQuadFace(sIndex, masterMap, false, forceOp);
        }
        else
        {
            coupledMesh& recvMesh = recvMeshes_[pI];

            // Bisect the slave face
            slaveMap =
            (
                recvMesh.subMesh().bisectQuadFace
                (
                    sIndex,
                    masterMap,
                    false,
                    forceOp
                )
            );
        }

        // Turn coupledModification back on.
        setCoupledModification();

        // The final operation has to succeed.
        if (slaveMap.type() <= 0)
        {
            FatalErrorIn
            (
                "\n"
                "const changeMap dynamicTopoFvMesh::bisectQuadFace\n"
                "(\n"
                "    const label fIndex,\n"
                "    const changeMap& masterMap,\n"
                "    bool checkOnly,\n"
                "    bool forceOp"
                ")\n"
            )
                << "Coupled topo-change for slave failed."
                << " Master face: " << fIndex << nl
                << " Slave face: " << sIndex << nl
                << " Patch index: " << pI << nl
                << " Type: " << slaveMap.type() << nl
                << abort(FatalError);
        }

        // Save index and patch for posterity
        slaveMap.index() = sIndex;
        slaveMap.patchIndex() = pI;
    }

    // Are we performing only checks?
    if (checkOnly)
    {
        map.type() = 1;
        return map;
    }

    if (debug > 1)
    {
        Pout<< nl << nl
            << "Face: " << fIndex
            << ": " << faces_[fIndex]
            << " is to be bisected. " << endl;

        if (debug > 2)
        {
            Pout<< " On SubMesh: " << isSubMesh_ << nl;
            Pout<< " coupledModification: " << coupledModification_ << nl;

            const polyBoundaryMesh& boundary = boundaryMesh();

            label epIndex = whichPatch(fIndex);

            Pout<< " Patch: ";

            if (epIndex == -1)
            {
                Pout<< "Internal" << endl;
            }
            else
            if (epIndex < boundary.size())
            {
                Pout<< boundary[epIndex].name() << endl;
            }
            else
            {
                Pout<< " New patch: " << epIndex << endl;
            }

            Pout<< "Cell[0]: " << c0 << ": " << oldCells[0] << endl;

            forAll(oldCells[0], faceI)
            {
                const labelList& fE = faceEdges_[oldCells[0][faceI]];

                Pout<< oldCells[0][faceI] << ": "
                    << faces_[oldCells[0][faceI]]
                    << " fE: " << fE
                    << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Pout<< '\t' << fE[edgeI]
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

    // Add the points to the map. Since this might require master mapping,
    // first check to see if a slave is being bisected.
    if (masterMap.addedPointList().size())
    {
        const List<objectMap>& pMap =
        (
            masterMap.addedPointList()
        );

        edge firstEdge
        (
            pMap[0].masterObjects()[0],
            pMap[0].masterObjects()[1]
        );

        edge secondEdge
        (
            pMap[1].masterObjects()[0],
            pMap[1].masterObjects()[1]
        );

        if (firstEdge == commonEdges[0] && secondEdge == commonEdges[1])
        {
            // Update in conventional order
            map.addPoint(newPointIndex[0]);
            map.addPoint(newPointIndex[1]);
        }
        else
        if (firstEdge == commonEdges[1] && secondEdge == commonEdges[0])
        {
            // Update in reverse order
            map.addPoint(newPointIndex[1]);
            map.addPoint(newPointIndex[0]);
        }
        else
        {
            // We have a problem
            FatalErrorIn
            (
                "\n"
                "const changeMap "
                "dynamicTopoFvMesh::bisectQuadFace\n"
                "(\n"
                "    const label fIndex,\n"
                "    const changeMap& masterMap,\n"
                "    bool checkOnly,\n"
                "    bool forceOp\n"
                ")\n"
            )
                << "Coupled topo-change for slave failed."
                << " firstEdge: " << firstEdge << nl
                << " secondEdge: " << secondEdge << nl
                << " commonEdges[0]: " << commonEdges[0] << nl
                << " commonEdges[1]: " << commonEdges[1] << nl
                << abort(FatalError);
        }
    }
    else
    {
        map.addPoint(newPointIndex[0]);
        map.addPoint(newPointIndex[1]);
    }

    // Add a new prism cell to the end of the list.
    // Currently invalid, but will be updated later.
    newCellIndex[0] = insertCell(newCells[0], lengthScale_[c0]);

    // Add this cell to the map.
    map.addCell(newCellIndex[0]);

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
        Pout<< "Modified face: " << fIndex
            << ": " << faces_[fIndex] << endl;

        if (debug > 2)
        {
            Pout<< "Common edges: " << nl
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
            newCellIndex[0],
            tmpQFEdges
        )
    );

    // Add this face to the map.
    map.addFace(newFaceIndex[0]);

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
            -1,
            tmpTFEdges
        )
    );

    // Add this face to the map.
    map.addFace(newFaceIndex[1]);

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
            -1,
            tmpTFEdges
        )
    );

    // Add this face to the map.
    map.addFace(newFaceIndex[2]);

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

    // Add this edge to the map.
    map.addEdge(newEdgeIndex[1]);

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

    // Add this edge to the map.
    map.addEdge(newEdgeIndex[2]);

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
                -1,
                tmpQFEdges
            )
        );

        // Add this face to the map.
        map.addFace(newFaceIndex[3]);

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

        // Add this edge to the map.
        map.addEdge(newEdgeIndex[0]);

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

        // Add this edge to the map.
        map.addEdge(newEdgeIndex[3]);

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

        // Add this edge to the map.
        map.addEdge(newEdgeIndex[4]);

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
            Pout<< "Modified Cell[0]: "
                << c0 << ": " << oldCells[0] << endl;

            forAll(oldCells[0], faceI)
            {
                const labelList& fE = faceEdges_[oldCells[0][faceI]];

                Pout<< oldCells[0][faceI]
                    << ": " << faces_[oldCells[0][faceI]]
                    << " fE: " << fE
                    << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Pout<< '\t' << fE[edgeI]
                        << ": " << edges_[fE[edgeI]]
                        << " eF: " << eF
                        << endl;
                }
            }

            Pout<< "New Cell[0]: " << newCellIndex[0]
                << ": " << newCells[0] << endl;

            forAll(newCells[0], faceI)
            {
                const labelList& fE = faceEdges_[newCells[0][faceI]];

                Pout<< newCells[0][faceI]
                    << ": " << faces_[newCells[0][faceI]]
                    << " fE: " << fE
                    << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Pout<< '\t' << fE[edgeI]
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

        // Add a new prism cell to the end of the list.
        // Currently invalid, but will be updated later.
        newCellIndex[1] = insertCell(newCells[1], lengthScale_[c1]);

        // Add this cell to the map.
        map.addCell(newCellIndex[1]);

        if (debug > 2)
        {
            Pout<< "Cell[1]: " << c1 << ": " << oldCells[1] << endl;

            forAll(oldCells[1], faceI)
            {
                const labelList& fE = faceEdges_[oldCells[1][faceI]];

                Pout<< oldCells[1][faceI] << ": "
                    << faces_[oldCells[1][faceI]]
                    << " fE: " << fE
                    << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Pout<< '\t' << fE[edgeI]
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
                newCellIndex[1],
                tmpQFEdges
            )
        );

        // Add this face to the map.
        map.addFace(newFaceIndex[3]);

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
            Pout<< "Common edges: " << nl
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
                newCellIndex[1],
                tmpQFEdges
            )
        );

        // Add this face to the map.
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
                -1,
                tmpTFEdges
            )
        );

        // Add this face to the map.
        map.addFace(newFaceIndex[5]);

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
                -1,
                tmpTFEdges
            )
        );

        // Add this face to the map.
        map.addFace(newFaceIndex[6]);

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

        // Add this edge to the map.
        map.addEdge(newEdgeIndex[0]);

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

        // Add this edge to the map.
        map.addEdge(newEdgeIndex[3]);

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

        // Add this edge to the map.
        map.addEdge(newEdgeIndex[4]);

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

        // Add this edge to the map.
        map.addEdge(newEdgeIndex[5]);

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

        // Add this edge to the map.
        map.addEdge(newEdgeIndex[6]);

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
            Pout<< nl << "Modified Cell[0]: "
                << c0 << ": " << oldCells[0] << endl;

            forAll(oldCells[0], faceI)
            {
                const labelList& fE = faceEdges_[oldCells[0][faceI]];

                Pout<< oldCells[0][faceI]
                    << ": " << faces_[oldCells[0][faceI]]
                    << " fE: " << fE
                    << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Pout<< '\t' << fE[edgeI]
                        << ": " << edges_[fE[edgeI]]
                        << " eF: " << eF
                        << endl;
                }
            }

            Pout<< "New Cell[0]: "
                << newCellIndex[0] << ": " << newCells[0] << endl;

            forAll(newCells[0], faceI)
            {
                const labelList& fE = faceEdges_[newCells[0][faceI]];

                Pout<< newCells[0][faceI] << ": "
                    << faces_[newCells[0][faceI]]
                    << " fE: " << fE
                    << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Pout<< '\t' << fE[edgeI]
                        << ": " << edges_[fE[edgeI]]
                        << " eF: " << eF
                        << endl;
                }
            }

            Pout<< nl << "Modified Cell[1]: "
                << c1 << ": " << oldCells[1] << endl;

            forAll(oldCells[1], faceI)
            {
                const labelList& fE = faceEdges_[oldCells[1][faceI]];

                Pout<< oldCells[1][faceI] << ": "
                    << faces_[oldCells[1][faceI]]
                    << " fE: " << fE
                    << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Pout<< '\t' << fE[edgeI]
                        << ": " << edges_[fE[edgeI]]
                        << " eF: " << eF
                        << endl;
                }
            }

            Pout<< "New Cell[1]: "
                << newCellIndex[1] << ": " << newCells[1] << endl;

            forAll(newCells[1], faceI)
            {
                const labelList& fE = faceEdges_[newCells[1][faceI]];

                Pout<< newCells[1][faceI] << ": "
                    << faces_[newCells[1][faceI]]
                    << " fE: " << fE
                    << endl;

                forAll(fE, edgeI)
                {
                    const labelList& eF = edgeFaces_[fE[edgeI]];

                    Pout<< '\t' << fE[edgeI]
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

    if (coupledModification_)
    {
        // Alias for convenience...
        const changeMap& slaveMap = slaveMaps[0];

        const label pI = slaveMap.patchIndex();

        // Fetch the appropriate coupleMap
        const coupleMap* cMapPtr = NULL;

        if (localCouple && !procCouple)
        {
            cMapPtr = &(patchCoupling_[pI].map());
        }
        else
        if (procCouple && !localCouple)
        {
            cMapPtr = &(recvMeshes_[pI].map());
        }

        // Alias for convenience
        const coupleMap& cMap = *cMapPtr;

        // Add the new points to the coupling map
        const List<objectMap>& apList = slaveMap.addedPointList();

        if (bisectingSlave)
        {
            // Update reverse pointMap
            cMap.mapMaster
            (
                coupleMap::POINT,
                newPointIndex[0],
                apList[0].index()
            );

            cMap.mapMaster
            (
                coupleMap::POINT,
                newPointIndex[1],
                apList[1].index()
            );

            // Update pointMap
            cMap.mapSlave
            (
                coupleMap::POINT,
                apList[0].index(),
                newPointIndex[0]
            );

            cMap.mapSlave
            (
                coupleMap::POINT,
                apList[1].index(),
                newPointIndex[1]
            );
        }
        else
        {
            // Update pointMap
            cMap.mapSlave
            (
                coupleMap::POINT,
                newPointIndex[0],
                apList[0].index()
            );

            cMap.mapSlave
            (
                coupleMap::POINT,
                newPointIndex[1],
                apList[1].index()
            );

            // Update reverse pointMap
            cMap.mapMaster
            (
                coupleMap::POINT,
                apList[0].index(),
                newPointIndex[0]
            );

            cMap.mapMaster
            (
                coupleMap::POINT,
                apList[1].index(),
                newPointIndex[1]
            );
        }

        // Create a master/slave entry for the new face on the patch.
        FixedList<bool, 2> foundMatch(false);
        FixedList<label, 2> checkFaces(-1);

        // Fill in the faces to check for...
        checkFaces[0] = fIndex;
        checkFaces[1] = newFaceIndex[3];

        const List<objectMap>& afList = slaveMap.addedFaceList();

        forAll (checkFaces, indexI)
        {
            const face& mFace = faces_[checkFaces[indexI]];

            label sFaceIndex = -1;

            forAll(afList, sfI)
            {
                const face* facePtr = NULL;

                if (localCouple && !procCouple)
                {
                    facePtr = &(faces_[afList[sfI].index()]);
                }
                else
                if (procCouple && !localCouple)
                {
                    const dynamicTopoFvMesh& mesh = recvMeshes_[pI].subMesh();

                    facePtr = &(mesh.faces_[afList[sfI].index()]);
                }

                const face& tFace = *facePtr;
                FixedList<bool, 4> cP(false);

                forAll(mFace, pointI)
                {
                    const Map<label>& pointMap =
                    (
                        cMap.entityMap(coupleMap::POINT)
                    );

                    if (tFace.which(pointMap[mFace[pointI]]) > -1)
                    {
                        cP[pointI] = true;
                    }
                }

                if (cP[0] && cP[1] && cP[2] && cP[3])
                {
                    sFaceIndex = afList[sfI].index();
                    foundMatch[indexI] = true;
                    break;
                }
            }

            if (foundMatch[indexI])
            {
                if (bisectingSlave)
                {
                    cMap.mapMaster
                    (
                        coupleMap::FACE,
                        checkFaces[indexI],
                        sFaceIndex
                    );

                    cMap.mapSlave
                    (
                        coupleMap::FACE,
                        sFaceIndex,
                        checkFaces[indexI]
                    );
                }
                else
                {
                    cMap.mapSlave
                    (
                        coupleMap::FACE,
                        checkFaces[indexI],
                        sFaceIndex
                    );

                    cMap.mapMaster
                    (
                        coupleMap::FACE,
                        sFaceIndex,
                        checkFaces[indexI]
                    );
                }
            }
            else
            {
                FatalErrorIn
                (
                    "\n"
                    "const changeMap "
                    "dynamicTopoFvMesh::bisectQuadFace\n"
                    "(\n"
                    "    const label fIndex,\n"
                    "    const changeMap& masterMap,\n"
                    "    bool checkOnly,\n"
                    "    bool forceOp\n"
                    ")\n"
                )
                    << "Failed to build coupled maps." << nl
                    << " masterFace: " << mFace << nl
                    << " Added slave faces: " << afList << nl
                    << abort(FatalError);
            }
        }

        // Push operation into coupleMap
        cMap.pushOperation
        (
            slaveMap.index(),
            coupleMap::BISECTION
        );
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
    status(TOTAL_BISECTIONS)++;

    // Increment surface-counter
    if (c1 == -1)
    {
        // Do not update stats for processor patches
        if (!processorCoupledEntity(fIndex))
        {
            status(SURFACE_BISECTIONS)++;
        }
    }

    // Increment the number of modifications
    status(TOTAL_MODIFICATIONS)++;

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
//    -3: Bisection failed since edge was on a noRefinement patch.
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
    //      Update faceEdges and edgeFaces information

    // For 2D meshes, perform face-bisection
    if (is2D())
    {
        return bisectQuadFace(eIndex, changeMap(), checkOnly);
    }

    // Figure out which thread this is...
    label tIndex = self();

    // Prepare the changeMaps
    changeMap map;
    List<changeMap> slaveMaps;
    bool bisectingSlave = false;

    if
    (
        (status(TOTAL_MODIFICATIONS) > maxModifications_) &&
        (maxModifications_ > -1)
    )
    {
        // Reached the max allowable topo-changes.
        stack(tIndex).clear();

        return map;
    }

    // Check if edgeRefinements are to be avoided on patch.
    const labelList& eF = edgeFaces_[eIndex];

    forAll(eF, fI)
    {
        label fPatch = whichPatch(eF[fI]);

        if (baseMesh_.lengthEstimator().checkRefinementPatch(fPatch))
        {
            map.type() = -3;

            return map;
        }
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
            << abort(FatalError);
    }

    // If coupled modification is set, and this is a
    // master edge, bisect its slaves as well.
    bool localCouple = false, procCouple = false;

    if (coupledModification_)
    {
        const edge& eCheck = edges_[eIndex];

        const label edgeEnum = coupleMap::EDGE;

        // Is this a locally coupled edge (either master or slave)?
        if (locallyCoupledEntity(eIndex, true))
        {
            localCouple = true;
            procCouple = false;
        }
        else
        if (processorCoupledEntity(eIndex))
        {
            procCouple = true;
            localCouple = false;
        }

        if (localCouple && !procCouple)
        {
            // Determine the slave index.
            label sIndex = -1, pIndex = -1;

            forAll(patchCoupling_, patchI)
            {
                if (patchCoupling_(patchI))
                {
                    const coupleMap& cMap = patchCoupling_[patchI].map();

                    if ((sIndex = cMap.findSlave(edgeEnum, eIndex)) > -1)
                    {
                        pIndex = patchI;

                        break;
                    }

                    // The following bit happens only during the sliver
                    // exudation process, since slave edges are
                    // usually not added to the coupled edge-stack.
                    if ((sIndex = cMap.findMaster(edgeEnum, eIndex)) > -1)
                    {
                        pIndex = patchI;

                        // Notice that we are bisecting a slave edge first.
                        bisectingSlave = true;

                        break;
                    }
                }
            }

            if (sIndex == -1)
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
                    << "Coupled maps were improperly specified." << nl
                    << " Slave index not found for: " << nl
                    << " Edge: " << eIndex << ": " << edges_[eIndex] << nl
                    << abort(FatalError);
            }
            else
            {
                // If we've found the slave, size up the list
                meshOps::sizeUpList
                (
                    changeMap(),
                    slaveMaps
                );

                // Save index and patch for posterity
                slaveMaps[0].index() = sIndex;
                slaveMaps[0].patchIndex() = pIndex;
            }

            if (debug > 1)
            {
                Pout<< nl << " >> Bisecting slave edge: " << sIndex
                    << " for master edge: " << eIndex << endl;
            }
        }
        else
        if (procCouple && !localCouple)
        {
            // If this is a new entity, bail out for now.
            // This will be handled at the next time-step.
            if (eIndex >= nOldEdges_)
            {
                return map;
            }

            // Check slaves
            forAll(procIndices_, pI)
            {
                // Fetch reference to subMeshes
                const coupledMesh& sendMesh = sendMeshes_[pI];
                const coupledMesh& recvMesh = recvMeshes_[pI];

                const coupleMap& scMap = sendMesh.map();
                const coupleMap& rcMap = recvMesh.map();

                // If this edge was sent to a lower-ranked
                // processor, skip it.
                if
                (
                    priority
                    (
                        procIndices_[pI],
                        lessOp<label>(),
                        Pstream::myProcNo()
                    )
                )
                {
                    if (scMap.reverseEntityMap(edgeEnum).found(eIndex))
                    {
                        if (debug > 3)
                        {
                            Pout<< "Edge: " << eIndex
                                << "::" << eCheck
                                << " was sent to proc: "
                                << procIndices_[pI]
                                << ", so bailing out."
                                << endl;
                        }

                        return map;
                    }
                }

                label sIndex = -1;

                if ((sIndex = rcMap.findSlave(edgeEnum, eIndex)) > -1)
                {
                    if (debug > 3)
                    {
                        Pout<< "Checking slave edge: " << sIndex
                            << " on proc: " << procIndices_[pI]
                            << " for master edge: " << eIndex
                            << endl;
                    }

                    // Check if a lower-ranked processor is
                    // handling this edge
                    if
                    (
                        priority
                        (
                            procIndices_[pI],
                            lessOp<label>(),
                            Pstream::myProcNo()
                        )
                    )
                    {
                        if (debug > 3)
                        {
                            Pout<< "Edge: " << eIndex
                                << " is handled by proc: "
                                << procIndices_[pI]
                                << ", so bailing out."
                                << endl;
                        }

                        return map;
                    }

                    label curIndex = slaveMaps.size();

                    // Size up the list
                    meshOps::sizeUpList
                    (
                        changeMap(),
                        slaveMaps
                    );

                    // Save index and patch for posterity
                    slaveMaps[curIndex].index() = sIndex;
                    slaveMaps[curIndex].patchIndex() = pI;
                }
            }
        }
        else
        {
            // Something's wrong with coupling maps
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
                << "Coupled maps were improperly specified." << nl
                << " localCouple: " << localCouple << nl
                << " procCouple: " << procCouple << nl
                << " Edge: " << eIndex << ": " << edges_[eIndex] << nl
                << abort(FatalError);
        }

        // Temporarily turn off coupledModification
        unsetCoupledModification();

        // First check the master for bisection feasibility.
        changeMap masterMap = bisectEdge(eIndex, true, forceOp);

        // Turn it back on.
        setCoupledModification();

        // Master couldn't perform bisection
        if (masterMap.type() <= 0)
        {
            return masterMap;
        }

        // Now check each of the slaves for bisection feasibility
        forAll(slaveMaps, slaveI)
        {
            // Alias for convenience...
            changeMap& slaveMap = slaveMaps[slaveI];

            label sIndex = slaveMap.index();
            label pI = slaveMap.patchIndex();

            // If the edge is mapped onto itself, skip check
            // (occurs for cyclic edges)
            if ((sIndex == eIndex) && localCouple)
            {
                continue;
            }

            // Temporarily turn off coupledModification
            unsetCoupledModification();

            // Test the slave edge
            if (localCouple)
            {
                slaveMap = bisectEdge(sIndex, true, forceOp);
            }
            else
            if (procCouple)
            {
                dynamicTopoFvMesh& sMesh = recvMeshes_[pI].subMesh();

                if (debug > 3)
                {
                    Pout<< "Checking slave edge: " << sIndex
                        << "::" << sMesh.edges_[sIndex]
                        << " on proc: " << procIndices_[pI]
                        << " for master edge: " << eIndex
                        << endl;
                }

                slaveMap = sMesh.bisectEdge(sIndex, true, forceOp);
            }

            // Turn it back on.
            setCoupledModification();

            if (slaveMap.type() <= 0)
            {
                // Slave couldn't perform bisection.
                map.type() = -2;

                return map;
            }

            // Save index and patch for posterity
            slaveMap.index() = sIndex;
            slaveMap.patchIndex() = pI;
        }

        // Next bisect each slave edge..
        forAll(slaveMaps, slaveI)
        {
            // Alias for convenience...
            changeMap& slaveMap = slaveMaps[slaveI];

            label sIndex = slaveMap.index();
            label pI = slaveMap.patchIndex();

            // If the edge is mapped onto itself, skip modification
            // (occurs for cyclic edges)
            if ((sIndex == eIndex) && localCouple)
            {
                continue;
            }

            // Temporarily turn off coupledModification
            unsetCoupledModification();

            // Bisect the slave edge
            if (localCouple)
            {
                slaveMap = bisectEdge(sIndex, false, forceOp);
            }
            else
            {
                dynamicTopoFvMesh& sMesh = recvMeshes_[pI].subMesh();

                slaveMap = sMesh.bisectEdge(sIndex, false, forceOp);
            }

            // Turn it back on.
            setCoupledModification();

            // The final operation has to succeed.
            if (slaveMap.type() <= 0)
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
                    << "Coupled topo-change for slave failed." << nl
                    << " Edge: " << eIndex << ": " << edges_[eIndex] << nl
                    << " Slave index: " << sIndex << nl
                    << " Patch index: " << pI << nl
                    << " Type: " << slaveMap.type() << nl
                    << abort(FatalError);
            }

            // Save index and patch for posterity
            slaveMap.index() = sIndex;
            slaveMap.patchIndex() = pI;
        }
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
        if (debug > 3 && isSubMesh_)
        {
            Pout<< "  Slave edge: " << eIndex
                << " can be bisected."
                << endl;
        }

        map.type() = 1;
        return map;
    }

    // Update number of surface bisections, if necessary.
    if (whichEdgePatch(eIndex) > -1)
    {
        // Do not update stats for processor patches
        if (!processorCoupledEntity(eIndex))
        {
            status(SURFACE_BISECTIONS)++;
        }
    }

    // Hull variables
    face tmpTriFace(3);
    edge origEdge(edges_[eIndex]);
    labelList tmpEdgeFaces(3,-1);
    labelList tmpIntEdgeFaces(4,-1);
    labelList tmpFaceEdges(3,-1);

    // Build vertexHull for this edge
    labelList vertexHull;
    buildVertexHull(eIndex, vertexHull);

    // Size up the hull lists
    label m = vertexHull.size();
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
        vertexHull,
        edgeHull,
        faceHull,
        cellHull,
        ringEntities
    );

    if (debug > 1)
    {
        Pout<< nl << nl
            << "Edge: " << eIndex
            << ": " << origEdge
            << " is to be bisected. " << endl;

        Pout<< " On SubMesh: " << isSubMesh_ << nl;
        Pout<< " coupledModification: " << coupledModification_ << nl;

        label epIndex = whichEdgePatch(eIndex);

        const polyBoundaryMesh& boundary = boundaryMesh();

        Pout<< " Patch: ";

        if (epIndex == -1)
        {
            Pout<< "Internal" << endl;
        }
        else
        if (epIndex < boundary.size())
        {
            Pout<< boundary[epIndex].name() << endl;
        }
        else
        {
            Pout<< " New patch: " << epIndex << endl;
        }

        // Write out VTK files prior to change
        //  - Using old-points for convenience in post-processing
        if (debug > 3)
        {
            writeVTK
            (
                Foam::name(eIndex)
              + '(' + Foam::name(origEdge[0])
              + ',' + Foam::name(origEdge[1]) + ')'
              + "_Bisect_0",
                cellHull,
                3, false, true
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

    if (locallyCoupledEntity(eIndex, true))
    {
        nePatch = locallyCoupledEdgePatch(eIndex);
    }
    else
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
            labelList(faceHull.size(),-1)
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
                insertCell
                (
                    newCell,
                    lengthScale_[cellHull[indexI]]
                )
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
                    addedCellIndices[indexI],
                    tmpFaceEdges
                )
            );

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
                        -1,
                        labelList(3, -1)
                    )
                );

                // Add this face to the map.
                map.addFace(addedFaceIndices[indexI]);

                // Configure edgeFaces
                tmpEdgeFaces[0] = faceHull[indexI];
                tmpEdgeFaces[1] = addedIntFaceIndices[indexI];
                tmpEdgeFaces[2] = addedFaceIndices[indexI];

                // Add an edge
                addedEdgeIndices[indexI] =
                (
                    insertEdge
                    (
                        whichPatch(faceHull[indexI]),
                        edge(newPointIndex,vertexHull[indexI]),
                        tmpEdgeFaces
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

                // Add the faceEdges entry
                faceEdges_[addedFaceIndices[indexI]] = tmpFaceEdges;

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
                        addedCellIndices[indexI],
                        labelList(3, -1)
                    )
                );

                // Add this face to the map.
                map.addFace(addedFaceIndices[indexI]);

                // Configure edgeFaces
                tmpIntEdgeFaces[0] = faceHull[indexI];
                tmpIntEdgeFaces[1] = addedIntFaceIndices[indexI];
                tmpIntEdgeFaces[2] = addedFaceIndices[indexI];
                tmpIntEdgeFaces[3] = addedIntFaceIndices[prevI];

                // Add an internal edge
                addedEdgeIndices[indexI] =
                (
                    insertEdge
                    (
                        -1,
                        edge(newPointIndex,vertexHull[indexI]),
                        tmpIntEdgeFaces
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

                // Add the faceEdges entry
                faceEdges_[addedFaceIndices[indexI]] = tmpFaceEdges;

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
                        addedCellIndices[indexI],
                        labelList(3, -1)
                    )
                );

                // Add this face to the map.
                map.addFace(addedFaceIndices[0]);

                // Configure edgeFaces
                tmpIntEdgeFaces[0] = faceHull[0];
                tmpIntEdgeFaces[1] = addedIntFaceIndices[0];
                tmpIntEdgeFaces[2] = addedFaceIndices[0];
                tmpIntEdgeFaces[3] = addedIntFaceIndices[indexI];

                // Add an internal edge
                addedEdgeIndices[0] =
                (
                    insertEdge
                    (
                        -1,
                        edge(newPointIndex,vertexHull[0]),
                        tmpIntEdgeFaces
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

                // Add the faceEdges entry
                faceEdges_[addedFaceIndices[0]] = tmpFaceEdges;

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
                    -1,
                    labelList(3, -1)
                )
            );

            // Add this face to the map.
            map.addFace(addedFaceIndices[indexI]);

            // Configure edgeFaces
            tmpEdgeFaces[0] = addedFaceIndices[indexI];
            tmpEdgeFaces[1] = addedIntFaceIndices[prevI];
            tmpEdgeFaces[2] = faceHull[indexI];

            // Add an edge
            addedEdgeIndices[indexI] =
            (
                insertEdge
                (
                    whichPatch(faceHull[indexI]),
                    edge(newPointIndex,vertexHull[indexI]),
                    tmpEdgeFaces
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

            // Add the faceEdges entry
            faceEdges_[addedFaceIndices[indexI]] = tmpFaceEdges;

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

    // If modification is coupled, update mapping info.
    if (coupledModification_)
    {
        // Build a list of boundary edges / faces for mapping
        dynamicLabelList checkEdges(8), checkFaces(4);

        const labelList& oeFaces = edgeFaces_[eIndex];
        const labelList& neFaces = edgeFaces_[newEdgeIndex];

        forAll(oeFaces, faceI)
        {
            FixedList<label, 2> check(-1);

            check[0] = oeFaces[faceI];
            check[1] = neFaces[faceI];

            forAll(check, indexI)
            {
                label cPatch = whichPatch(check[indexI]);

                if (localCouple && !procCouple)
                {
                    if (!locallyCoupledEntity(check[indexI], true, false, true))
                    {
                        continue;
                    }

                    // Check for cyclics
                    if (boundaryMesh()[cPatch].type() == "cyclic")
                    {
                        // Check if this is a master face
                        const coupleMap& cM = patchCoupling_[cPatch].map();
                        const Map<label>& fM = cM.entityMap(coupleMap::FACE);

                        // Only add master faces
                        // (belonging to the first half)
                        //  - Only check[0] will be present,
                        //    so check for that alone
                        if (!fM.found(check[0]))
                        {
                            continue;
                        }
                    }
                }
                else
                if (procCouple && !localCouple)
                {
                    if (getNeighbourProcessor(cPatch) == -1)
                    {
                        continue;
                    }
                }

                // Add face and its edges for checking
                if (findIndex(checkFaces, check[indexI]) == -1)
                {
                    // Add this face
                    checkFaces.append(check[indexI]);

                    const labelList& fEdges = faceEdges_[check[indexI]];

                    forAll(fEdges, edgeI)
                    {
                        if (findIndex(checkEdges, fEdges[edgeI]) == -1)
                        {
                            checkEdges.append(fEdges[edgeI]);
                        }
                    }
                }
            }
        }

        // Prepare a checklist
        boolList matchedFaces(checkFaces.size(), false);
        boolList matchedEdges(checkEdges.size(), false);

        // Output check entities
        if (debug > 4)
        {
            writeVTK
            (
                "checkEdges_" + Foam::name(eIndex),
                checkEdges, 1, false, true
            );

            writeVTK
            (
                "checkFaces_" + Foam::name(eIndex),
                checkFaces, 2, false, true
            );
        }

        forAll(slaveMaps, slaveI)
        {
            // Alias for convenience...
            const changeMap& slaveMap = slaveMaps[slaveI];

            const label pI = slaveMap.patchIndex();

            // Fetch the appropriate coupleMap / mesh
            const coupleMap* cMapPtr = NULL;
            const dynamicTopoFvMesh* sMeshPtr = NULL;

            if (localCouple && !procCouple)
            {
                cMapPtr = &(patchCoupling_[pI].map());
                sMeshPtr = this;
            }
            else
            if (procCouple && !localCouple)
            {
                cMapPtr = &(recvMeshes_[pI].map());
                sMeshPtr = &(recvMeshes_[pI].subMesh());
            }

            // Alias for convenience
            const coupleMap& cMap = *cMapPtr;
            const dynamicTopoFvMesh& sMesh = *sMeshPtr;

            // Add the new points to the coupling map
            const List<objectMap>& apList = slaveMap.addedPointList();

            // Fetch the slave point
            label slavePoint = -1;

            if ((slaveMap.index() == eIndex) && localCouple)
            {
                // Cyclic edge at axis
                slavePoint = newPointIndex;
            }
            else
            {
                slavePoint = apList[0].index();
            }

            // Map points
            if (bisectingSlave)
            {
                // Update reverse pointMap
                cMap.mapMaster
                (
                    coupleMap::POINT,
                    newPointIndex,
                    slavePoint
                );

                // Update pointMap
                cMap.mapSlave
                (
                    coupleMap::POINT,
                    slavePoint,
                    newPointIndex
                );
            }
            else
            {
                // Update pointMap
                cMap.mapSlave
                (
                    coupleMap::POINT,
                    newPointIndex,
                    slavePoint
                );

                // Update reverse pointMap
                cMap.mapMaster
                (
                    coupleMap::POINT,
                    slavePoint,
                    newPointIndex
                );
            }

            if (debug > 2)
            {
                Pout<< " Adding point: " << slavePoint
                    << " for point: " << newPointIndex
                    << endl;
            }

            // Obtain point maps
            const Map<label>& pointMap = cMap.entityMap(coupleMap::POINT);

            // Update face mapping
            const label faceEnum = coupleMap::FACE;

            // Obtain references
            Map<label>& faceMap = cMap.entityMap(faceEnum);
            Map<label>& rFaceMap = cMap.reverseEntityMap(faceEnum);

            forAll(checkFaces, faceI)
            {
                label mfIndex = checkFaces[faceI];

                const face& mF = faces_[mfIndex];

                label mfPatch = whichPatch(mfIndex);

                // Check for processor match
                label neiProc = getNeighbourProcessor(mfPatch);

                if (procCouple && !localCouple)
                {
                    if (neiProc != procIndices_[pI])
                    {
                        continue;
                    }
                }

                triFace cF
                (
                    pointMap.found(mF[0]) ? pointMap[mF[0]] : -1,
                    pointMap.found(mF[1]) ? pointMap[mF[1]] : -1,
                    pointMap.found(mF[2]) ? pointMap[mF[2]] : -1
                );

                // Skip mapping if all points were not found
                if (cF[0] == -1 || cF[1] == -1 || cF[2] == -1)
                {
                    continue;
                }

                // Fetch edges connected to the slave point
                const labelList& spEdges = sMesh.pointEdges_[slavePoint];

                forAll(spEdges, edgeI)
                {
                    label seIndex = spEdges[edgeI];

                    if (sMesh.whichEdgePatch(seIndex) == -1)
                    {
                        continue;
                    }

                    const labelList& seFaces = sMesh.edgeFaces_[seIndex];

                    forAll(seFaces, faceJ)
                    {
                        label sfIndex = seFaces[faceJ];

                        if (sMesh.whichPatch(sfIndex) == -1)
                        {
                            continue;
                        }

                        const face& sF = sMesh.faces_[sfIndex];

                        if
                        (
                            triFace::compare
                            (
                                triFace(sF[0], sF[1], sF[2]), cF
                            )
                        )
                        {
                            if (debug > 2)
                            {
                                word pN(boundaryMesh()[mfPatch].name());

                                Pout<< " Found face: " << sfIndex
                                    << " :: " << sF
                                    << " with mfIndex: " << mfIndex
                                    << " :: " << mF
                                    << " Patch: " << pN
                                    << endl;
                            }

                            if (rFaceMap.found(sfIndex))
                            {
                                rFaceMap[sfIndex] = mfIndex;
                            }
                            else
                            {
                                rFaceMap.insert(sfIndex, mfIndex);
                            }

                            if (faceMap.found(mfIndex))
                            {
                                faceMap[mfIndex] = sfIndex;
                            }
                            else
                            {
                                faceMap.insert(mfIndex, sfIndex);
                            }

                            matchedFaces[faceI] = true;

                            break;
                        }
                    }

                    if (matchedFaces[faceI])
                    {
                        break;
                    }
                }

                if ((debug > 4) && !matchedFaces[faceI])
                {
                    sMesh.writeVTK
                    (
                        "failedFacePoints_"
                      + Foam::name(mfIndex),
                        labelList(cF), 0, false, true
                    );

                    writeVTK
                    (
                        "checkFaces_" + Foam::name(eIndex),
                        checkFaces, 2, false, true
                    );

                    Pout<< " Failed to match face: "
                        << mfIndex << " :: " << mF
                        << " masterPatch: " << mfPatch
                        << " using comparison face: " << cF
                        << " on proc: " << procIndices_[pI]
                        << endl;
                }
            }

            // Update edge mapping
            const label edgeEnum = coupleMap::EDGE;

            // Obtain references
            Map<label>& edgeMap = cMap.entityMap(edgeEnum);
            Map<label>& rEdgeMap = cMap.reverseEntityMap(edgeEnum);

            forAll(checkEdges, edgeI)
            {
                label meIndex = checkEdges[edgeI];

                const edge& mE = edges_[meIndex];

                label mePatch = whichEdgePatch(meIndex);
                label neiProc = getNeighbourProcessor(mePatch);

                edge cE
                (
                    pointMap.found(mE[0]) ? pointMap[mE[0]] : -1,
                    pointMap.found(mE[1]) ? pointMap[mE[1]] : -1
                );

                // Skip mapping if all points were not found
                if (cE[0] == -1 || cE[1] == -1)
                {
                    continue;
                }

                // Fetch edges connected to the slave point
                const labelList& spEdges = sMesh.pointEdges_[cE[0]];

                forAll(spEdges, edgeJ)
                {
                    label seIndex = spEdges[edgeJ];

                    const edge& sE = sMesh.edges_[seIndex];

                    if (sE == cE)
                    {
                        if (debug > 2)
                        {
                            Pout<< " Found edge: " << seIndex
                                << " :: " << sE
                                << " with meIndex: " << meIndex
                                << " :: " << mE
                                << " on proc: " << procIndices_[pI]
                                << endl;
                        }

                        // Update reverse map
                        if (rEdgeMap.found(seIndex))
                        {
                            rEdgeMap[seIndex] = meIndex;
                        }
                        else
                        {
                            rEdgeMap.insert(seIndex, meIndex);
                        }

                        // Update map
                        if (edgeMap.found(meIndex))
                        {
                            edgeMap[meIndex] = seIndex;
                        }
                        else
                        {
                            edgeMap.insert(meIndex, seIndex);
                        }

                        matchedEdges[edgeI] = true;

                        break;
                    }
                }

                if (!matchedEdges[edgeI])
                {
                    if (procCouple && !localCouple)
                    {
                        // Rare occassion where both points
                        // of the edge lie on processor, but
                        // not the edge itself.
                        if (neiProc != procIndices_[pI])
                        {
                            if (debug > 2)
                            {
                                Pout<< " Edge: " << meIndex
                                    << " :: " << mE
                                    << " with comparison: " << cE
                                    << " has points on processor: " << neiProc
                                    << " but no edge. Marking as matched."
                                    << endl;
                            }

                            matchedEdges[edgeI] = true;
                        }
                    }
                }

                if ((debug > 4) && !matchedEdges[edgeI])
                {
                    sMesh.writeVTK
                    (
                        "failedEdge_"
                      + Foam::name(meIndex),
                        labelList(cE), 0, false, true
                    );

                    writeVTK
                    (
                        "checkEdges_" + Foam::name(eIndex),
                        checkEdges, 1, false, true
                    );

                    Pout<< " Failed to match edge: "
                        << meIndex << " :: " << mE
                        << " using comparison edge: " << cE
                        << " on proc: " << procIndices_[pI]
                        << endl;
                }
            }

            // Push operation into coupleMap
            if (procCouple && !localCouple)
            {
                cMap.pushOperation
                (
                    slaveMap.index(),
                    coupleMap::BISECTION
                );
            }
        }

        // Ensure that all entities were matched
        label nFailFace = 0, nFailEdge = 0;

        forAll(matchedFaces, faceI)
        {
            if (!matchedFaces[faceI])
            {
                ++nFailFace;
            }
        }

        forAll(matchedEdges, edgeI)
        {
            if (!matchedEdges[edgeI])
            {
                ++nFailEdge;
            }
        }

        if (nFailFace || nFailEdge)
        {
            Pout<< " Failed to match all entities. " << nl
                << "  Faces: " << nFailFace << nl
                << "  Edges: " << nFailEdge << nl
                << abort(FatalError);
        }
    }

    if (debug > 2)
    {
        label bPatch = whichEdgePatch(eIndex);

        const polyBoundaryMesh& boundary = boundaryMesh();

        if (bPatch == -1)
        {
            Pout<< "Patch: Internal" << nl;
        }
        else
        if (bPatch < boundary.size())
        {
            Pout<< "Patch: " << boundary[bPatch].name() << nl;
        }
        else
        {
            Pout<< " New patch: " << bPatch << endl;
        }

        Pout<< "vertexHull: " << vertexHull << nl
            << "Edges: " << edgeHull << nl
            << "Faces: " << faceHull << nl
            << "Cells: " << cellHull << nl;

        Pout<< "Modified cells: " << nl;

        forAll(cellHull, cellI)
        {
            if (cellHull[cellI] == -1)
            {
                continue;
            }

            Pout<< cellHull[cellI] << ":: "
                << cells_[cellHull[cellI]]
                << nl;
        }

        Pout<< nl << "Added cells: " << nl;

        forAll(addedCellIndices, cellI)
        {
            if (addedCellIndices[cellI] == -1)
            {
                continue;
            }

            Pout<< addedCellIndices[cellI] << ":: "
                << cells_[addedCellIndices[cellI]] << nl
                << "lengthScale: " << lengthScale_[addedCellIndices[cellI]]
                << nl;
        }

        Pout<< nl << "Modified faces: " << nl;

        forAll(faceHull, faceI)
        {
            Pout<< faceHull[faceI] << ":: "
                << faces_[faceHull[faceI]] << ": "
                << owner_[faceHull[faceI]] << ": "
                << neighbour_[faceHull[faceI]] << " "
                << "faceEdges:: " << faceEdges_[faceHull[faceI]]
                << nl;
        }

        Pout<< nl << "Added faces: " << nl;

        forAll(addedFaceIndices, faceI)
        {
            Pout<< addedFaceIndices[faceI] << ":: "
                << faces_[addedFaceIndices[faceI]] << ": "
                << owner_[addedFaceIndices[faceI]] << ": "
                << neighbour_[addedFaceIndices[faceI]] << " "
                << "faceEdges:: " << faceEdges_[addedFaceIndices[faceI]]
                << nl;
        }

        forAll(addedIntFaceIndices, faceI)
        {
            if (addedIntFaceIndices[faceI] == -1)
            {
                continue;
            }

            Pout<< addedIntFaceIndices[faceI] << ":: "
                << faces_[addedIntFaceIndices[faceI]] << ": "
                << owner_[addedIntFaceIndices[faceI]] << ": "
                << neighbour_[addedIntFaceIndices[faceI]] << " "
                << "faceEdges:: " << faceEdges_[addedIntFaceIndices[faceI]]
                << nl;
        }

        Pout<< nl << "New edge:: " << newEdgeIndex
            << ": " << edges_[newEdgeIndex] << nl
            << " edgeFaces:: " << edgeFaces_[newEdgeIndex]
            << nl;

        Pout<< nl << "Added edges: " << nl;

        forAll(addedEdgeIndices, edgeI)
        {
            Pout<< addedEdgeIndices[edgeI]
                << ":: " << edges_[addedEdgeIndices[edgeI]] << nl
                << " edgeFaces:: " << edgeFaces_[addedEdgeIndices[edgeI]]
                << nl;
        }

        Pout<< "New Point:: " << newPointIndex << nl
            << "pointEdges:: " << pointEdges_[newPointIndex] << nl;

        // Flush buffer
        Pout<< endl;

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

            // Write out VTK files after change
            //  - Using old-points for convenience in post-processing
            writeVTK
            (
                Foam::name(eIndex)
              + '(' + Foam::name(origEdge[0])
              + ',' + Foam::name(origEdge[1]) + ')'
              + "_Bisect_1",
                newHull,
                3, false, true
            );
        }
    }

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    status(TOTAL_BISECTIONS)++;

    // Increment the number of modifications
    status(TOTAL_MODIFICATIONS)++;

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

    // Obtain a reference to this edge
    const edge& edgeToCheck = edges_[eIndex];

    // Obtain point references
    const point& a = points_[edgeToCheck[0]];
    const point& c = points_[edgeToCheck[1]];

    const point& aOld = oldPoints_[edgeToCheck[0]];
    const point& cOld = oldPoints_[edgeToCheck[1]];

    // Compute the mid-point of the edge
    point midPoint = 0.5 * (a + c);
    point oldPoint = 0.5 * (aOld + cOld);

    dynamicLabelList eCells(10);

    const labelList& eFaces = edgeFaces_[eIndex];

    // Accumulate cells connected to this edge
    forAll(eFaces, faceI)
    {
        label own = owner_[eFaces[faceI]];
        label nei = neighbour_[eFaces[faceI]];

        if (findIndex(eCells, own) == -1)
        {
            eCells.append(own);
        }

        if (nei == -1)
        {
            continue;
        }

        if (findIndex(eCells, nei) == -1)
        {
            eCells.append(nei);
        }
    }

    // Loop through all cells and compute quality
    forAll(eCells, cellI)
    {
        label cellIndex = eCells[cellI];

        const cell& checkCell = cells_[cellIndex];

        // Find two faces that don't contain the edge
        forAll(checkCell, faceI)
        {
            const face& checkFace = faces_[checkCell[faceI]];

            bool check0 = (findIndex(checkFace, edgeToCheck[0]) == -1);
            bool check1 = (findIndex(checkFace, edgeToCheck[1]) == -1);

            if ((check0 && !check1) || (!check0 && check1))
            {
                // Check orientation
                if (owner_[checkCell[faceI]] == cellIndex)
                {
                    cQuality =
                    (
                        tetMetric_
                        (
                            points_[checkFace[2]],
                            points_[checkFace[1]],
                            points_[checkFace[0]],
                            midPoint
                        )
                    );

                    oldVolume =
                    (
                        tetPointRef
                        (
                            oldPoints_[checkFace[2]],
                            oldPoints_[checkFace[1]],
                            oldPoints_[checkFace[0]],
                            oldPoint
                        ).mag()
                    );
                }
                else
                {
                    cQuality =
                    (
                        tetMetric_
                        (
                            points_[checkFace[0]],
                            points_[checkFace[1]],
                            points_[checkFace[2]],
                            midPoint
                        )
                    );

                    oldVolume =
                    (
                        tetPointRef
                        (
                            oldPoints_[checkFace[0]],
                            oldPoints_[checkFace[1]],
                            oldPoints_[checkFace[2]],
                            oldPoint
                        ).mag()
                    );
                }

                // Check if the quality is worse
                minQuality = Foam::min(cQuality, minQuality);
                minVolume = Foam::min(oldVolume, minVolume);
            }
        }
    }

    // Ensure that the mesh is valid
    if (minQuality < sliverThreshold_)
    {
        if (debug > 3 && minQuality < 0.0)
        {
            writeVTK(Foam::name(eIndex) + "_iCells", eCells);
        }

        if (debug > 2)
        {
            InfoIn
            (
                "scalar dynamicTopoFvMesh::computeBisectionQuality"
                "(const label eIndex) const"
            )
                << " Bisecting edge will fall below the"
                << " sliver threshold of: " << sliverThreshold_ << nl
                << " Edge: " << eIndex << ": " << edgeToCheck << nl
                << " Minimum Quality: " << minQuality << nl
                << " Minimum Volume: " << minVolume << nl
                << " Mid point: " << midPoint << nl
                << " Old point: " << oldPoint << nl
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


// Slice the mesh at a particular location
void dynamicTopoFvMesh::sliceMesh
(
    const labelPair& pointPair
)
{
    if (debug > 1)
    {
        Pout<< nl << nl
            << "Pair: " << pointPair
            << " is to be used for mesh slicing. " << endl;
    }

    label patchIndex = -1;
    scalar dx = 0.0;
    vector gCentre = vector::zero;
    FixedList<vector, 2> fC(vector::zero);

    if (is2D())
    {
        patchIndex = whichPatch(pointPair.first());

        // Fetch face centres
        fC[0] = faces_[pointPair.first()].centre(points_);
        fC[1] = faces_[pointPair.second()].centre(points_);
    }
    else
    {
        // Find the patch that the edge-vertex is connected to.
        const labelList& pEdges = pointEdges_[pointPair.first()];

        forAll(pEdges, edgeI)
        {
            if ((patchIndex = whichEdgePatch(pEdges[edgeI])) > -1)
            {
                break;
            }
        }

        fC[0] = points_[pointPair.first()];
        fC[1] = points_[pointPair.second()];
    }

    linePointRef lpr(fC[0], fC[1]);

    // Specify the centre.
    gCentre = lpr.centre();

    // Specify a search distance
    dx = lpr.mag();

    // Is this edge in the vicinity of a previous slice-point?
    if (lengthEstimator().checkOldSlices(gCentre))
    {
        if (debug > 1)
        {
            Pout<< nl << nl
                << "Pair: " << pointPair
                << " is too close to another slice point. "
                << endl;
        }

        // Too close to another slice-point. Bail out.
        return;
    }

    // Choose a box around the centre and scan all
    // surface entities that fall into this region.
    boundBox bBox
    (
        gCentre - vector(dx, dx, dx),
        gCentre + vector(dx, dx, dx)
    );

    vector p = vector::zero, N = vector::zero;
    Map<vector> checkPoints, surfFaces;
    Map<edge> checkEdges;

    if (is2D())
    {
        // Assign plane point / normal
        p = gCentre;

        vector gNorm = faces_[pointPair.first()].normal(points_);

        gNorm /= (mag(gNorm) + VSMALL);

        // Since this is 2D, assume XY-plane here.
        N = (gNorm ^ vector(0.0, 0.0, 1.0));
    }
    else
    {
        // Prepare surface points / edges for Dijkstra's algorithm
        for (label edgeI = nOldInternalEdges_; edgeI < nEdges_; edgeI++)
        {
            if (edgeFaces_[edgeI].empty())
            {
                continue;
            }

            if (whichEdgePatch(edgeI) == patchIndex)
            {
                const edge& surfaceEdge = edges_[edgeI];

                if
                (
                    (bBox.contains(points_[surfaceEdge[0]])) &&
                    (bBox.contains(points_[surfaceEdge[1]]))
                )
                {
                    checkEdges.insert(edgeI, surfaceEdge);

                    if (!checkPoints.found(surfaceEdge[0]))
                    {
                        checkPoints.insert
                        (
                            surfaceEdge[0],
                            points_[surfaceEdge[0]]
                        );
                    }

                    if (!checkPoints.found(surfaceEdge[1]))
                    {
                        checkPoints.insert
                        (
                            surfaceEdge[1],
                            points_[surfaceEdge[1]]
                        );
                    }

                    // Add surface faces as well.
                    const labelList& eFaces = edgeFaces_[edgeI];

                    forAll(eFaces, faceI)
                    {
                        if
                        (
                            (neighbour_[eFaces[faceI]] == -1) &&
                            (!surfFaces.found(eFaces[faceI]))
                        )
                        {
                            vector surfNorm =
                            (
                                faces_[eFaces[faceI]].normal(points_)
                            );

                            surfFaces.insert(eFaces[faceI], surfNorm);
                        }
                    }
                }
            }
        }

        if (debug > 1)
        {
            Pout<< nl << nl
                << " Point [0]: " << points_[pointPair.first()] << nl
                << " Point [1]: " << points_[pointPair.second()] << endl;

            if (debug > 3)
            {
                writeVTK("slicePoints", checkPoints.toc(), 0);
                writeVTK("sliceEdges", checkEdges.toc(), 1);
            }
        }

        // Find the shortest path using Dijkstra's algorithm.
        Map<label> shortestPath;

        bool foundPath =
        (
            meshOps::Dijkstra
            (
                checkPoints,
                checkEdges,
                pointPair.first(),
                pointPair.second(),
                shortestPath
            )
        );

        // First normalize all face-normals
        forAllIter(Map<vector>, surfFaces, sIter)
        {
            sIter() /= (mag(sIter()) + VSMALL);
        }

        if (foundPath)
        {
            // Next, take cross-products with every other
            // vector in the list, and accumulate.
            forAllIter(Map<vector>, surfFaces, sIterI)
            {
                forAllIter(Map<vector>, surfFaces, sIterJ)
                {
                    if (sIterI.key() != sIterJ.key())
                    {
                        vector n = (sIterI() ^ sIterJ());

                        n /= (mag(n) + VSMALL);

                        // Reverse the vector if necessary
                        if ((N & n) < 0.0)
                        {
                            n = -n;
                        }

                        N += n;
                    }
                }
            }

            N /= (mag(N) + VSMALL);

            // Obtain point position
            p = gCentre;
        }
        else
        {
            // Probably a membrane-type configuration.
            labelHashSet checkCells;

            // Prepare a bounding cylinder with radius dx.
            forAllIter(Map<vector>, surfFaces, sIter)
            {
                const face& thisFace = faces_[sIter.key()];

                if (thisFace.which(pointPair.first()) > -1)
                {
                    N += sIter();
                }
            }

            // Normalize and reverse.
            N /= -(mag(N) + VSMALL);

            vector a0 = points_[pointPair.first()];
            vector a1 = points_[pointPair.second()];
            scalar dist = mag(a1 - a0);

            forAll(cells_, cellI)
            {
                if (cells_[cellI].empty())
                {
                    continue;
                }

                scalar v = 0.0;
                vector x = vector::zero;

                meshOps::cellCentreAndVolume
                (
                    cellI,
                    points_,
                    faces_,
                    cells_,
                    owner_,
                    x,
                    v
                );

                vector rx = (x - a0);
                vector ra = (rx & N)*N;

                // Check if point falls off cylinder ends.
                if (mag(ra) > dist || mag(ra) < 0.0)
                {
                    continue;
                }

                vector r = (rx - ra);

                // Check if the magnitude of 'r' is within radius.
                if (mag(r) < dx)
                {
                    checkCells.insert(cellI);
                }
            }

            labelList cList = checkCells.toc();

            if (debug > 1)
            {
                Pout<< "Dijkstra's algorithm could not find a path." << endl;

                if (debug > 3)
                {
                    writeVTK("checkCells", cList, 3);
                }
            }

            changeMap sliceMap =
            (
                removeCells
                (
                    cList,
                    patchIndex,
                    "Slice_"
                  + Foam::name(pointPair.first()) + '_'
                  + Foam::name(pointPair.second())
                )
            );

            if (debug)
            {
                checkConnectivity(10);
            }

            // Accumulate a list of projection points
            labelHashSet pPoints;

            const List<objectMap>& afList = sliceMap.addedFaceList();

            forAll(afList, faceI)
            {
                const face& thisFace = faces_[afList[faceI].index()];

                forAll(thisFace, pointI)
                {
                    if (!pPoints.found(thisFace[pointI]))
                    {
                        pPoints.insert(thisFace[pointI]);
                    }
                }
            }

            // Now project points in a series of intermediate steps.

            // Add an entry to sliceBoxes.
            lengthEstimator().appendBox(bBox);

            return;
        }
    }

    if (debug > 1)
    {
        Pout<< nl << nl
            << " Plane point: " << p << nl
            << " Plane normal: " << N << endl;
    }

    // Mark cells and interior faces that fall
    // within the bounding box.
    labelHashSet checkCells, checkFaces, splitFaces;
    Map<bool> cellColors;

    forAll(faces_, faceI)
    {
        if (faces_[faceI].empty())
        {
            continue;
        }

        if (is2D() && faces_[faceI].size() == 3)
        {
            continue;
        }

        vector fCentre = faces_[faceI].centre(points_);

        FixedList<label, 2> cellsToCheck(-1);
        cellsToCheck[0] = owner_[faceI];
        cellsToCheck[1] = neighbour_[faceI];

        if (bBox.contains(fCentre) && cellsToCheck[1] != -1)
        {
            // Add this internal face to the list.
            checkFaces.insert(faceI);

            scalar volume = 0.0;
            vector centre = vector::zero;

            forAll(cellsToCheck, cellI)
            {
                if (!checkCells.found(cellsToCheck[cellI]))
                {
                    meshOps::cellCentreAndVolume
                    (
                        cellsToCheck[cellI],
                        points_,
                        faces_,
                        cells_,
                        owner_,
                        centre,
                        volume
                    );

                    checkCells.insert(cellsToCheck[cellI]);

                    if (((centre - p) & N) > 0.0)
                    {
                        cellColors.insert(cellsToCheck[cellI], true);
                    }
                    else
                    {
                        cellColors.insert(cellsToCheck[cellI], false);
                    }
                }
            }
        }
    }

    // Prepare a list of internal faces for mesh splitting.
    forAllIter(labelHashSet, checkFaces, fIter)
    {
        if
        (
            cellColors[owner_[fIter.key()]]
         != cellColors[neighbour_[fIter.key()]]
        )
        {
            splitFaces.insert(fIter.key());
        }

        // Loop through all points (and associated pointEdges)
        // for this face, and check if connected cells are also
        // present in the checkCells/cellColors list
        if (is2D())
        {
            const labelList& fEdges = faceEdges_[fIter.key()];

            forAll(fEdges, edgeI)
            {
                const labelList& eFaces = edgeFaces_[fEdges[edgeI]];

                forAll(eFaces, faceI)
                {
                    label own = owner_[eFaces[faceI]];
                    label nei = neighbour_[eFaces[faceI]];

                    if (!checkCells.found(own))
                    {
                        scalar volume = 0.0;
                        vector centre = vector::zero;

                        meshOps::cellCentreAndVolume
                        (
                            own,
                            points_,
                            faces_,
                            cells_,
                            owner_,
                            centre,
                            volume
                        );

                        checkCells.insert(own);

                        if (((centre - p) & N) > 0.0)
                        {
                            cellColors.insert(own, true);
                        }
                        else
                        {
                            cellColors.insert(own, false);
                        }
                    }

                    if (!checkCells.found(nei) && nei != -1)
                    {
                        scalar volume = 0.0;
                        vector centre = vector::zero;

                        meshOps::cellCentreAndVolume
                        (
                            nei,
                            points_,
                            faces_,
                            cells_,
                            owner_,
                            centre,
                            volume
                        );

                        checkCells.insert(nei);

                        if (((centre - p) & N) > 0.0)
                        {
                            cellColors.insert(nei, true);
                        }
                        else
                        {
                            cellColors.insert(nei, false);
                        }
                    }
                }
            }
        }
        else
        {
            const face& faceToCheck = faces_[fIter.key()];

            forAll(faceToCheck, pointI)
            {
                const labelList& pEdges = pointEdges_[faceToCheck[pointI]];

                forAll(pEdges, edgeI)
                {
                    const labelList& eFaces = edgeFaces_[pEdges[edgeI]];

                    forAll(eFaces, faceI)
                    {
                        label own = owner_[eFaces[faceI]];
                        label nei = neighbour_[eFaces[faceI]];

                        if (!checkCells.found(own))
                        {
                            scalar volume = 0.0;
                            vector centre = vector::zero;

                            meshOps::cellCentreAndVolume
                            (
                                own,
                                points_,
                                faces_,
                                cells_,
                                owner_,
                                centre,
                                volume
                            );

                            checkCells.insert(own);

                            if (((centre - p) & N) > 0.0)
                            {
                                cellColors.insert(own, true);
                            }
                            else
                            {
                                cellColors.insert(own, false);
                            }
                        }

                        if (!checkCells.found(nei) && nei != -1)
                        {
                            scalar volume = 0.0;
                            vector centre = vector::zero;

                            meshOps::cellCentreAndVolume
                            (
                                nei,
                                points_,
                                faces_,
                                cells_,
                                owner_,
                                centre,
                                volume
                            );

                            checkCells.insert(nei);

                            if (((centre - p) & N) > 0.0)
                            {
                                cellColors.insert(nei, true);
                            }
                            else
                            {
                                cellColors.insert(nei, false);
                            }
                        }
                    }
                }
            }
        }
    }

    if (debug > 3)
    {
        writeVTK("splitFaces", splitFaces.toc(), 2);
        writeVTK("checkCells", checkCells.toc(), 3);
    }

    // Pass this info into the splitInternalFaces routine.
    splitInternalFaces
    (
        patchIndex,
        splitFaces.toc(),
        cellColors
    );

    // Add an entry to sliceBoxes.
    lengthEstimator().appendBox(bBox);
}


// Add cell layer above specified patch
const changeMap dynamicTopoFvMesh::addCellLayer
(
    const label patchID
)
{
    changeMap map;

    // Maps for added entities
    Map<label> addedPoints;
    Map<label> addedHEdges, addedVEdges, currentVEdges;
    Map<label> addedHFaces, addedVFaces, currentVFaces;
    Map<labelPair> addedCells;

    // Allocate a list of patch faces
    dynamicLabelList patchFaces(patchSizes_[patchID]);

    // Loop through all patch faces and create a cell for each
    for (label faceI = nOldInternalFaces_; faceI < faces_.size(); faceI++)
    {
        label pIndex = whichPatch(faceI);

        if (pIndex != patchID)
        {
            continue;
        }

        // Add face to the list
        patchFaces.append(faceI);

        // Add a new cell
        label cIndex = owner_[faceI];
        scalar newLengthScale = -1.0;
        const cell& ownCell = cells_[cIndex];

        if (edgeRefinement_)
        {
            newLengthScale = lengthScale_[cIndex];
        }

        label newCellIndex =
        (
            insertCell
            (
                cell(ownCell.size()),
                newLengthScale
            )
        );

        // Update maps
        map.addCell(newCellIndex, labelList(1, cIndex));
        addedCells.insert(cIndex, labelPair(newCellIndex, 0));
    }

    labelList mP(2, -1);

    forAll(patchFaces, indexI)
    {
        label faceI = patchFaces[indexI];
        label cIndex = owner_[faceI];

        // Fetch appropriate face / cell
        //  - Make copies, since holding references
        //    to data within this loop is unsafe.
        const face bFace = faces_[faceI];
        const cell ownCell = cells_[cIndex];

        // Configure a new face for insertion
        face newHFace(bFace);
        labelList newHFaceEdges(bFace.size(), -1);

        // Get the opposing face from the cell
        oppositeFace oFace = ownCell.opposingFace(faceI, faces_);

        if (!oFace.found())
        {
            // Something's wrong here.
            FatalErrorIn
            (
                "const changeMap dynamicTopoFvMesh::addCellLayer"
                "(const label patchID)"
            )
                << " Face: " << faceI << " :: " << bFace << nl
                << " has no opposing face in cell: "
                << cIndex << " :: " << ownCell << nl
                << abort(FatalError);
        }

        // Create points
        forAll(bFace, pointI)
        {
            label pIndex = bFace[pointI];

            // Skip if we've added this already
            if (addedPoints.found(pIndex))
            {
                continue;
            }

            // Set master points
            mP[0] = pIndex;
            mP[1] = oFace[pointI];

            label newPointIndex =
            (
                insertPoint
                (
                    0.5 * (points_[mP[0]] + points_[mP[1]]),
                    oldPoints_[mP[0]],
                    mP
                )
            );

            // Update maps
            map.addPoint(newPointIndex, mP);
            addedPoints.insert(pIndex, newPointIndex);
        }

        // Fetch faceEdges from opposite faces.
        //  - Make copies, since holding references is unsafe
        const labelList bfEdges = faceEdges_[faceI];
        const labelList ofEdges = faceEdges_[oFace.oppositeIndex()];

        // Create edges for each edge of the new horizontal face
        forAll(bfEdges, edgeI)
        {
            label beIndex = bfEdges[edgeI];

            // Skip if we've added this already
            if (addedHEdges.found(beIndex))
            {
                // Update face edges for the new horizontal face
                newHFaceEdges[edgeI] = addedHEdges[beIndex];

                continue;
            }

            // Configure the new edge
            label oeIndex = -1;
            const edge bEdge = edges_[beIndex];

            // Build an edge for comparison
            edge cEdge
            (
                oFace[bFace.which(bEdge[0])],
                oFace[bFace.which(bEdge[1])]
            );

            forAll(ofEdges, edgeJ)
            {
                const edge& ofEdge = edges_[ofEdges[edgeJ]];

                if (cEdge == ofEdge)
                {
                    oeIndex = ofEdges[edgeJ];
                    break;
                }
            }

            if (oeIndex < 0)
            {
                FatalErrorIn
                (
                    "const changeMap dynamicTopoFvMesh::addCellLayer"
                    "(const label patchID)"
                )
                    << " Could not find comparison edge: " << cEdge
                    << " for edge: " << bEdge
                    << abort(FatalError);
            }

            // Fetch patch information
            label hEdgePatch = whichEdgePatch(oeIndex);

            // Set indices
            edge newHEdge
            (
                addedPoints[bEdge[0]],
                addedPoints[bEdge[1]]
            );

            // Insert a new edge with empty edgeFaces
            label newHEdgeIndex =
            (
                insertEdge
                (
                    hEdgePatch,
                    newHEdge,
                    labelList(0)
                )
            );

            // Update maps
            map.addEdge(newHEdgeIndex);
            addedHEdges.insert(beIndex, newHEdgeIndex);

            // Update face edges for the new horizontal face
            newHFaceEdges[edgeI] = newHEdgeIndex;

            // Add a new vertical face for this edge
            label vFaceIndex = -1;

            // Find a vertical face that contains both edges
            const labelList& beFaces = edgeFaces_[beIndex];

            forAll(beFaces, faceJ)
            {
                const labelList& testEdges = faceEdges_[beFaces[faceJ]];

                if
                (
                    (findIndex(testEdges, beIndex) > -1) &&
                    (findIndex(testEdges, oeIndex) > -1)
                )
                {
                    vFaceIndex = beFaces[faceJ];
                    break;
                }
            }

            if (vFaceIndex < 0)
            {
                FatalErrorIn
                (
                    "const changeMap dynamicTopoFvMesh::addCellLayer"
                    "(const label patchID)"
                )
                    << " Could not find an appropriate vertical face"
                    << " containing edges: " << oeIndex
                    << " and " << beIndex
                    << abort(FatalError);
            }

            // Find two vertical edges on this face
            const labelList& vfEdges = faceEdges_[vFaceIndex];

            forAll(vfEdges, edgeJ)
            {
                const edge& vfEdge = edges_[vfEdges[edgeJ]];

                forAll(bEdge, i)
                {
                    if (vfEdge == edge(bEdge[i], cEdge[i]))
                    {
                        // Skip if we've added this already
                        if (addedVEdges.found(bEdge[i]))
                        {
                            continue;
                        }

                        label veIndex = vfEdges[edgeJ];

                        // Fetch edge patch information
                        label vEdgePatch = whichEdgePatch(veIndex);

                        // Set indices
                        edge newVEdge
                        (
                            bEdge[i],
                            addedPoints[bEdge[i]]
                        );

                        // Insert a new edge with empty edgeFaces
                        label newVEdgeIndex =
                        (
                            insertEdge
                            (
                                vEdgePatch,
                                newVEdge,
                                labelList(0)
                            )
                        );

                        // Update maps
                        map.addEdge(newVEdgeIndex);
                        addedVEdges.insert(bEdge[i], newVEdgeIndex);

                        // Note edge indices for later renumbering
                        currentVEdges.insert(bEdge[i], veIndex);
                    }
                }
            }

            // Configure the new vertical face
            face newVFace(faces_[vFaceIndex]);
            label newOwner = -1, newNeighbour = -1;

            label oldOwner = owner_[vFaceIndex];
            label oldNeighbour = neighbour_[vFaceIndex];

            // Fetch owner / neighbour
            newOwner = addedCells[oldOwner].first();

            if (oldNeighbour > -1)
            {
                newNeighbour = addedCells[oldNeighbour].first();
            }

            // Replace point indices on the new face
            forAll(bEdge, i)
            {
                meshOps::replaceLabel
                (
                    cEdge[i],
                    addedPoints[bEdge[i]],
                    newVFace
                );
            }

            // Note face indices for later renumbering
            currentVFaces.insert(beIndex, vFaceIndex);

            // Check if reversal is necessary
            if ((newNeighbour < newOwner) && (newNeighbour > -1))
            {
                // Flip face
                newVFace = newVFace.reverseFace();

                // Swap addressing
                Foam::Swap(newOwner, newNeighbour);
                Foam::Swap(oldOwner, oldNeighbour);
            }

            // Configure faceEdges for the new vertical face
            labelList newVFaceEdges(4, -1);

            newVFaceEdges[0] = beIndex;
            newVFaceEdges[1] = newHEdgeIndex;
            newVFaceEdges[2] = addedVEdges[bEdge[0]];
            newVFaceEdges[3] = addedVEdges[bEdge[1]];

            // Add the new vertical face
            label newVFaceIndex =
            (
                insertFace
                (
                    whichPatch(vFaceIndex),
                    newVFace,
                    newOwner,
                    newNeighbour,
                    newVFaceEdges
                )
            );

            // Update maps
            map.addFace(newVFaceIndex, labelList(1, vFaceIndex));
            addedVFaces.insert(beIndex, newVFaceIndex);

            // Update face count on the new cells
            cells_[newOwner][addedCells[oldOwner].second()++] =
            (
                newVFaceIndex
            );

            if (newNeighbour > -1)
            {
                cells_[newNeighbour][addedCells[oldNeighbour].second()++] =
                (
                    newVFaceIndex
                );
            }

            // Size up edgeFaces for each edge
            forAll(newVFaceEdges, edgeJ)
            {
                label vfeIndex = newVFaceEdges[edgeJ];

                meshOps::sizeUpList
                (
                    newVFaceIndex,
                    edgeFaces_[vfeIndex]
                );
            }
        }

        // Add a new interior face, with identical orientation
        forAll(newHFace, pointI)
        {
            newHFace[pointI] = addedPoints[newHFace[pointI]];
        }

        // Add the new horizontal face
        label newHFaceIndex =
        (
            insertFace
            (
                -1,
                newHFace,
                cIndex,
                addedCells[cIndex].first(),
                newHFaceEdges
            )
        );

        // Update maps
        map.addFace(newHFaceIndex, labelList(1, faceI));
        addedHFaces.insert(faceI, newHFaceIndex);

        // Replace index on the old cell
        meshOps::replaceLabel
        (
            faceI,
            newHFaceIndex,
            cells_[cIndex]
        );

        // Update face count on the new cell
        label newCellIndex = addedCells[cIndex].first();

        // Modify owner for the existing boundary face
        owner_[faceI] = newCellIndex;

        cells_[newCellIndex][addedCells[cIndex].second()++] = faceI;
        cells_[newCellIndex][addedCells[cIndex].second()++] = newHFaceIndex;

        // Size up edgeFaces for each edge
        forAll(newHFaceEdges, edgeI)
        {
            label hfeIndex = newHFaceEdges[edgeI];

            meshOps::sizeUpList
            (
                newHFaceIndex,
                edgeFaces_[hfeIndex]
            );
        }
    }

    // Renumber vertical edges
    forAllConstIter(Map<label>, currentVEdges, eIter)
    {
        // Fetch reference to edge
        edge& curEdge = edges_[eIter()];

        if (curEdge[0] == eIter.key())
        {
            curEdge[0] = addedPoints[eIter.key()];
        }

        if (curEdge[1] == eIter.key())
        {
            curEdge[1] = addedPoints[eIter.key()];
        }

        // Size down pointEdges
        if (is3D())
        {
            meshOps::sizeDownList
            (
                eIter(),
                pointEdges_[eIter.key()]
            );

            meshOps::sizeUpList
            (
                eIter(),
                pointEdges_[addedPoints[eIter.key()]]
            );
        }
    }

    // Renumber vertical faces
    forAllConstIter(Map<label>, currentVFaces, fIter)
    {
        // Fetch reference to existing edge
        const edge& bEdge = edges_[fIter.key()];

        // Replace point indices on vertical face
        forAll(bEdge, i)
        {
            meshOps::replaceLabel
            (
                bEdge[i],
                addedPoints[bEdge[i]],
                faces_[fIter()]
            );
        }

        // Replace edge on the existing vertical face
        meshOps::replaceLabel
        (
            fIter.key(),
            addedHEdges[fIter.key()],
            faceEdges_[fIter()]
        );

        // Remove old face on existing boundary edge
        meshOps::sizeDownList
        (
            fIter(),
            edgeFaces_[fIter.key()]
        );

        // Add old face to new horizontal edge
        meshOps::sizeUpList
        (
            fIter(),
            edgeFaces_[addedHEdges[fIter.key()]]
        );
    }

    // Now that all old / new cells possess correct connectivity,
    // compute mapping information.
    const List<objectMap>& afList = map.addedFaceList();

    forAll(afList, indexI)
    {
        label parent = afList[indexI].masterObjects()[0];

        if (whichPatch(afList[indexI].index()) == -1)
        {
            // Interior faces get default mapping
            if (whichPatch(parent) == -1)
            {
                setFaceMapping(parent);
            }

            setFaceMapping(afList[indexI].index());
        }
    }

    // Set the flag
    topoChangeFlag_ = true;

    // Specify that the operation was successful
    map.type() = 1;

    // Return the changeMap
    return map;
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

        if (is3D())
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

        if (is3D())
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
            // Make an identical faceEdges entry.
            // This will be renumbered once new edges are added.
            labelList newFaceEdges(faceEdges_[internalFaces[faceI]]);

            newFaceIndex[indexI] =
            (
                insertFace
                (
                    patchIndex,
                    newFace[indexI],
                    newOwner[indexI],
                    -1,
                    newFaceEdges
                )
            );

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
    // We'll deal with correcting edgeFaces later.
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
                newEdgeFaces
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
                    newEdgeFaces
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

    if (is2D())
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
