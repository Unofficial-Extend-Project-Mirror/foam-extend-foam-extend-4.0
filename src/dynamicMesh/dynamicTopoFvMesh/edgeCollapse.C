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

// Method to collapse a quad-face in 2D
// - Returns a changeMap with a type specifying:
//    -3: Collapse failed since face was on a noRefinement patch.
//    -1: Collapse failed since max number of topo-changes was reached.
//     0: Collapse could not be performed.
//     1: Collapsed to first node.
//     2: Collapsed to second node.
//     3: Collapse to mid-point.
// - overRideCase is used to force a certain collapse configuration.
//    -1: Use this value to let collapseQuadFace decide a case.
//     1: Force collapse to first node.
//     2: Force collapse to second node.
//     3: Force collapse to mid-point.
// - checkOnly performs a feasibility check and returns without modifications.
const changeMap dynamicTopoFvMesh::collapseQuadFace
(
    const label fIndex,
    label overRideCase,
    bool checkOnly,
    bool forceOp
)
{
    // Figure out which thread this is...
    label tIndex = self();

    // Prepare the changeMaps
    changeMap map;
    List<changeMap> slaveMaps;
    bool collapsingSlave = false;

    if
    (
        (status(TOTAL_MODIFICATIONS) > maxModifications_)
     && (maxModifications_ > -1)
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
            "const changeMap "
            "dynamicTopoFvMesh::collapseQuadFace\n"
            "(\n"
            "    const label fIndex,\n"
            "    label overRideCase,\n"
            "    bool checkOnly,\n"
            "    bool forceOp\n"
            ")\n"
        )
            << " Invalid index: " << fIndex
            << abort(FatalError);
    }

    // Define the edges on the face to be collapsed
    FixedList<edge,4> checkEdge(edge(-1, -1));
    FixedList<label,4> checkEdgeIndex(-1);

    // Define checkEdges
    getCheckEdges(fIndex, (*this), map, checkEdge, checkEdgeIndex);

    // Determine the common vertices for the first and second edges
    label cv0 = checkEdge[1].commonVertex(checkEdge[0]);
    label cv1 = checkEdge[1].commonVertex(checkEdge[3]);
    label cv2 = checkEdge[2].commonVertex(checkEdge[0]);
    label cv3 = checkEdge[2].commonVertex(checkEdge[3]);

    // If coupled modification is set, and this is a
    // master face, collapse its slaves first.
    bool localCouple = false, procCouple = false;

    if (coupledModification_)
    {
        const face& fCheck = faces_[fIndex];

        const label faceEnum = coupleMap::FACE;
        const label pointEnum = coupleMap::POINT;

        // Is this a locally coupled edge (either master or slave)?
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
                    // exudation process, since slave edges are
                    // usually not added to the coupled edge-stack.
                    if ((sIndex = cMap.findMaster(faceEnum, fIndex)) > -1)
                    {
                        pIndex = patchI;

                        // Notice that we are collapsing a slave edge first.
                        collapsingSlave = true;

                        break;
                    }
                }
            }

            if (sIndex == -1)
            {
                FatalErrorIn
                (
                    "\n"
                    "const changeMap "
                    "dynamicTopoFvMesh::collapseQuadFace\n"
                    "(\n"
                    "    const label fIndex,\n"
                    "    label overRideCase,\n"
                    "    bool checkOnly,\n"
                    "    bool forceOp\n"
                    ")\n"
                )
                    << "Coupled maps were improperly specified." << nl
                    << " Slave index not found for: " << nl
                    << " Face: " << fIndex << ": " << fCheck << nl
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
                Pout<< nl << " >> Collapsing slave face: " << sIndex
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
                // Fetch reference to subMeshes
                const coupledMesh& sendMesh = sendMeshes_[pI];
                const coupledMesh& recvMesh = recvMeshes_[pI];

                const coupleMap& scMap = sendMesh.map();
                const coupleMap& rcMap = recvMesh.map();

                // If this face was sent to a lower-ranked
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
                    if (scMap.reverseEntityMap(faceEnum).found(fIndex))
                    {
                        if (debug > 3)
                        {
                            Pout<< "Face: " << fIndex
                                << "::" << fCheck
                                << " was sent to proc: "
                                << procIndices_[pI]
                                << ", so bailing out."
                                << endl;
                        }

                        return map;
                    }
                }

                label sIndex = -1;

                if ((sIndex = rcMap.findSlave(faceEnum, fIndex)) > -1)
                {
                    // Check if a lower-ranked processor is
                    // handling this face
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
                                << "::" << fCheck
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
                else
                if
                (
                    (
                        rcMap.findSlave(pointEnum, cv0) > -1 &&
                        rcMap.findSlave(pointEnum, cv1) > -1
                    )
                 || (
                        rcMap.findSlave(pointEnum, cv2) > -1 &&
                        rcMap.findSlave(pointEnum, cv3) > -1
                    )
                )
                {
                    // An edge-only coupling exists.

                    // Check if a lower-ranked processor is
                    // handling this face
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
                             Pout<< "Face edge on: " << fIndex
                                 << "::" << fCheck
                                 << " is handled by proc: "
                                 << procIndices_[pI]
                                 << ", so bailing out."
                                 << endl;
                         }

                         return map;
                    }

                    label p0 = rcMap.findSlave(pointEnum, cv0);
                    label p1 = rcMap.findSlave(pointEnum, cv1);

                    label p2 = rcMap.findSlave(pointEnum, cv2);
                    label p3 = rcMap.findSlave(pointEnum, cv3);

                    edge cEdge(-1, -1);
                    label edgeCouple = -1;

                    if ((p0 > -1 && p1 > -1) && (p2 == -1 && p3 == -1))
                    {
                        cEdge[0] = p0;
                        cEdge[1] = p1;

                        edgeCouple = 1;
                        sIndex = readLabel(map.lookup("firstEdge"));
                    }
                    else
                    if ((p0 == -1 && p1 == -1) && (p2 > -1 && p3 > -1))
                    {
                        cEdge[0] = p2;
                        cEdge[1] = p3;

                        edgeCouple = 2;
                        sIndex = readLabel(map.lookup("secondEdge"));
                    }

                    label curIndex = slaveMaps.size();

                    // Size up the list
                    meshOps::sizeUpList
                    (
                        changeMap(),
                        slaveMaps
                    );

                    // Unfortunately, since no edge maps are
                    // available in 2D, loop through boundary
                    // faces on subMesh and get the right edge
                    label seIndex = -1;

                    const dynamicTopoFvMesh& sMesh = recvMesh.subMesh();

                    for
                    (
                        label faceI = sMesh.nOldInternalFaces_;
                        faceI < sMesh.faceEdges_.size();
                        faceI++
                    )
                    {
                        const labelList& fEdges = sMesh.faceEdges_[faceI];

                        forAll(fEdges, edgeI)
                        {
                            if (sMesh.edges_[fEdges[edgeI]] == cEdge)
                            {
                                seIndex = fEdges[edgeI];
                                break;
                            }
                        }

                        if (seIndex > -1)
                        {
                            break;
                        }
                    }

                    if (seIndex == -1)
                    {
                        Pout<< "Face edge on: " << fIndex
                            << "::" << fCheck
                            << " sIndex: " << sIndex
                            << " could not be found."
                            << abort(FatalError);
                    }

                    // Save index and patch for posterity
                    //  - Negate the index to signify edge coupling
                    slaveMaps[curIndex].index() = -seIndex;
                    slaveMaps[curIndex].patchIndex() = pI;

                    // Save edgeCouple as well, so that
                    // another map comparison is avoided.
                    slaveMaps[curIndex].type() = edgeCouple;
                }
            }
        }
        else
        {
            // Something's wrong with coupling maps
            FatalErrorIn
            (
                "\n"
                "const changeMap "
                "dynamicTopoFvMesh::collapseQuadFace\n"
                "(\n"
                "    const label fIndex,\n"
                "    label overRideCase,\n"
                "    bool checkOnly,\n"
                "    bool forceOp\n"
                ")\n"
            )
                << "Coupled maps were improperly specified." << nl
                << " localCouple: " << localCouple << nl
                << " procCouple: " << procCouple << nl
                << " Face: " << fIndex << ": " << fCheck << nl
                << abort(FatalError);
        }

        // Temporarily turn off coupledModification
        unsetCoupledModification();

        // Test the master face for collapse, and decide on a case
        changeMap masterMap = collapseQuadFace(fIndex, -1, true, forceOp);

        // Turn it back on.
        setCoupledModification();

        // Master couldn't perform collapse
        if (masterMap.type() <= 0)
        {
            return masterMap;
        }

        // For edge-only coupling, define the points for checking
        List<FixedList<point, 2> > slaveMoveNewPoint(slaveMaps.size());
        List<FixedList<point, 2> > slaveMoveOldPoint(slaveMaps.size());

        // Now check each of the slaves for collapse feasibility
        forAll(slaveMaps, slaveI)
        {
            // Alias for convenience...
            changeMap& slaveMap = slaveMaps[slaveI];

            label slaveOverRide = -1;
            label sIndex = slaveMap.index();
            label pI = slaveMap.patchIndex();
            const coupleMap* cMapPtr = NULL;

            if (localCouple)
            {
                cMapPtr = &(patchCoupling_[pI].map());
            }
            else
            if (procCouple)
            {
                const dynamicTopoFvMesh& sMesh =
                (
                    recvMeshes_[pI].subMesh()
                );

                cMapPtr = &(recvMeshes_[pI].map());

                if (debug > 3)
                {
                    if (sIndex < 0)
                    {
                        Pout<< "Checking slave edge: " << mag(sIndex)
                            << "::" << sMesh.edges_[mag(sIndex)]
                            << " on proc: " << procIndices_[pI]
                            << " for master face: " << fIndex
                            << " using collapseCase: " << masterMap.type()
                            << endl;
                    }
                    else
                    {
                        Pout<< "Checking slave face: " << sIndex
                            << "::" << sMesh.faces_[sIndex]
                            << " on proc: " << procIndices_[pI]
                            << " for master face: " << fIndex
                            << " using collapseCase: " << masterMap.type()
                            << endl;
                    }
                }
            }

            // Define an overRide case for the slave
            FixedList<edge, 2> mEdge(edge(-1, -1)), sEdge(edge(-1, -1));

            if (collapsingSlave)
            {
                const Map<label>& rPointMap =
                (
                    cMapPtr->reverseEntityMap(pointEnum)
                );

                if (sIndex < 0)
                {
                    if (slaveMap.type() == 1)
                    {
                        mEdge[0][0] = rPointMap[cv0];
                        mEdge[0][1] = rPointMap[cv1];
                    }
                    else
                    if (slaveMap.type() == 2)
                    {
                        mEdge[0][0] = rPointMap[cv2];
                        mEdge[0][1] = rPointMap[cv3];
                    }
                }
                else
                {
                    mEdge[0][0] = rPointMap[cv0];
                    mEdge[0][1] = rPointMap[cv1];

                    mEdge[1][0] = rPointMap[cv2];
                    mEdge[1][1] = rPointMap[cv3];
                }
            }
            else
            {
                const Map<label>& pointMap =
                (
                    cMapPtr->entityMap(pointEnum)
                );

                if (sIndex < 0)
                {
                    if (slaveMap.type() == 1)
                    {
                        mEdge[0][0] = pointMap[cv0];
                        mEdge[0][1] = pointMap[cv1];
                    }
                    else
                    if (slaveMap.type() == 2)
                    {
                        mEdge[0][0] = pointMap[cv2];
                        mEdge[0][1] = pointMap[cv3];
                    }
                }
                else
                {
                    mEdge[0][0] = pointMap[cv0];
                    mEdge[0][1] = pointMap[cv1];

                    mEdge[1][0] = pointMap[cv2];
                    mEdge[1][1] = pointMap[cv3];
                }
            }

            // Determine checkEdges for the slave
            FixedList<edge,4> slaveCheckEdge(edge(-1, -1));
            FixedList<label,4> slaveCheckEdgeIndex(-1);

            if (localCouple)
            {
                getCheckEdges
                (
                    sIndex,
                    (*this),
                    slaveMap,
                    slaveCheckEdge,
                    slaveCheckEdgeIndex
                );

                sEdge[0] = edges_[readLabel(slaveMap.lookup("firstEdge"))];
                sEdge[1] = edges_[readLabel(slaveMap.lookup("secondEdge"))];
            }
            else
            if (procCouple)
            {
                const dynamicTopoFvMesh& sMesh =
                (
                    recvMeshes_[pI].subMesh()
                );

                if (sIndex < 0)
                {
                    sEdge[0] = sMesh.edges_[mag(sIndex)];
                }
                else
                {
                    getCheckEdges
                    (
                        sIndex,
                        sMesh,
                        slaveMap,
                        slaveCheckEdge,
                        slaveCheckEdgeIndex
                    );

                    sEdge[0] =
                    (
                        sMesh.edges_[readLabel(slaveMap.lookup("firstEdge"))]
                    );

                    sEdge[1] =
                    (
                        sMesh.edges_[readLabel(slaveMap.lookup("secondEdge"))]
                    );
                }
            }

            // Compare edge orientations for edge-only coupling
            label compVal = -2;

            // Perform a topological comparison.
            switch (masterMap.type())
            {
                case 1:
                {
                    if (sIndex < 0)
                    {
                        slaveMoveNewPoint[slaveI][0] = points_[cv0];
                        slaveMoveNewPoint[slaveI][1] = points_[cv1];

                        slaveMoveOldPoint[slaveI][0] = oldPoints_[cv0];
                        slaveMoveOldPoint[slaveI][1] = oldPoints_[cv1];
                    }
                    else
                    if (mEdge[0] == sEdge[0])
                    {
                        slaveOverRide = 1;
                    }
                    else
                    if (mEdge[1] == sEdge[0])
                    {
                        slaveOverRide = 2;
                    }
                    else
                    {
                        // Write out for for post-processing
                        writeVTK("mFace_" + Foam::name(fIndex), fIndex, 2);

                        FatalErrorIn
                        (
                            "\n"
                            "const changeMap "
                            "dynamicTopoFvMesh::collapseQuadFace\n"
                            "(\n"
                            "    const label fIndex,\n"
                            "    label overRideCase,\n"
                            "    bool checkOnly,\n"
                            "    bool forceOp\n"
                            ")\n"
                        )
                            << "Coupled collapse failed." << nl
                            << "Masters: " << nl
                            << checkEdgeIndex[1] << ": "
                            << checkEdge[1] << nl
                            << checkEdgeIndex[2] << ": "
                            << checkEdge[2] << nl
                            << "Slaves: " << nl
                            << readLabel(slaveMap.lookup("firstEdge")) << ": "
                            << sEdge[0] << nl
                            << readLabel(slaveMap.lookup("secondEdge")) << ": "
                            << sEdge[1] << nl
                            << abort(FatalError);
                    }

                    break;
                }

                case 2:
                {
                    if (sIndex < 0)
                    {
                        slaveMoveNewPoint[slaveI][0] = points_[cv2];
                        slaveMoveNewPoint[slaveI][1] = points_[cv3];

                        slaveMoveOldPoint[slaveI][0] = oldPoints_[cv2];
                        slaveMoveOldPoint[slaveI][1] = oldPoints_[cv3];
                    }
                    else
                    if (mEdge[1] == sEdge[1])
                    {
                        slaveOverRide = 2;
                    }
                    else
                    if (mEdge[0] == sEdge[1])
                    {
                        slaveOverRide = 1;
                    }
                    else
                    {
                        // Write out for for post-processing
                        writeVTK("mFace_" + Foam::name(fIndex), fIndex, 2);

                        FatalErrorIn
                        (
                            "\n"
                            "const changeMap "
                            "dynamicTopoFvMesh::collapseQuadFace\n"
                            "(\n"
                            "    const label fIndex,\n"
                            "    label overRideCase,\n"
                            "    bool checkOnly,\n"
                            "    bool forceOp\n"
                            ")\n"
                        )
                            << "Coupled collapse failed." << nl
                            << "Masters: " << nl
                            << checkEdgeIndex[1] << ": "
                            << checkEdge[1] << nl
                            << checkEdgeIndex[2] << ": "
                            << checkEdge[2] << nl
                            << "Slaves: " << nl
                            << readLabel(slaveMap.lookup("firstEdge")) << ": "
                            << sEdge[0] << nl
                            << readLabel(slaveMap.lookup("secondEdge")) << ": "
                            << sEdge[1] << nl
                            << abort(FatalError);
                    }

                    break;
                }

                case 3:
                {
                    if (sIndex < 0)
                    {
                        slaveMoveNewPoint[slaveI][0] =
                        (
                            0.5 * (points_[cv0] + points_[cv2])
                        );

                        slaveMoveNewPoint[slaveI][1] =
                        (
                            0.5 * (points_[cv1] + points_[cv3])
                        );

                        slaveMoveOldPoint[slaveI][0] =
                        (
                            0.5 * (oldPoints_[cv0] + oldPoints_[cv2])
                        );

                        slaveMoveOldPoint[slaveI][1] =
                        (
                            0.5 * (oldPoints_[cv1] + oldPoints_[cv3])
                        );
                    }
                    else
                    {
                        overRideCase = 3;
                    }

                    break;
                }
            }

            if (sIndex < 0)
            {
                // Check edge orientation
                compVal = edge::compare(mEdge[0], sEdge[0]);

                // Swap components if necessary
                if (compVal == -1)
                {
                    Foam::Swap
                    (
                        slaveMoveNewPoint[slaveI][0],
                        slaveMoveNewPoint[slaveI][1]
                    );

                    Foam::Swap
                    (
                        slaveMoveOldPoint[slaveI][0],
                        slaveMoveOldPoint[slaveI][1]
                    );
                }
                else
                if (compVal != 1)
                {
                    FatalErrorIn
                    (
                        "\n"
                        "const changeMap "
                        "dynamicTopoFvMesh::collapseQuadFace\n"
                        "(\n"
                        "    const label fIndex,\n"
                        "    label overRideCase,\n"
                        "    bool checkOnly,\n"
                        "    bool forceOp\n"
                        ")\n"
                    )
                        << "Coupled topo-change for slave failed." << nl
                        << " Dissimilar edges: " << nl
                        << " mEdge: " << mEdge << nl
                        << " sEdge: " << sEdge << nl
                        << " masterMap type: " << masterMap.type() << nl
                        << abort(FatalError);
                }
            }

            // Temporarily turn off coupledModification
            unsetCoupledModification();

            // Test the slave face
            if (localCouple)
            {
                slaveMap =
                (
                    collapseQuadFace(sIndex, slaveOverRide, true, forceOp)
                );
            }
            else
            if (procCouple)
            {
                dynamicTopoFvMesh& sMesh = recvMeshes_[pI].subMesh();

                if (sIndex < 0)
                {
                    // Edge-based coupling

                    // Build a hull of cells and tri-faces
                    // that are connected to the slave edge
                    labelList slaveHullCells;
                    labelList slaveHullTriFaces;

                    meshOps::constructPrismHull
                    (
                        mag(sIndex),
                        sMesh.faces_,
                        sMesh.cells_,
                        sMesh.owner_,
                        sMesh.neighbour_,
                        sMesh.edgeFaces_,
                        slaveHullTriFaces,
                        slaveHullCells
                    );

                    bool infeasible = false;

                    // Keep track of resulting cell quality,
                    // if collapse is indeed feasible
                    scalar slaveCollapseQuality(GREAT);

                    // Check whether the collapse is possible.
                    if
                    (
                        checkCollapse
                        (
                            sMesh,
                            slaveHullTriFaces,
                            FixedList<label, 2>(-1),
                            FixedList<label, 2>(-1),
                            sEdge[0],
                            slaveMoveNewPoint[slaveI],
                            slaveMoveOldPoint[slaveI],
                            slaveCollapseQuality,
                            false,
                            forceOp
                        )
                    )
                    {
                        infeasible = true;
                    }

                    if (infeasible)
                    {
                        slaveMap.type() = 0;
                    }
                    else
                    {
                        slaveMap.type() = 1;
                    }
                }
                else
                {
                    // Edge-based coupling
                    slaveMap =
                    (
                        sMesh.collapseQuadFace
                        (
                            sIndex,
                            slaveOverRide,
                            true,
                            forceOp
                        )
                    );
                }
            }

            // Turn it back on.
            setCoupledModification();

            if (slaveMap.type() <= 0)
            {
                // Slave couldn't perform collapse.
                map.type() = -2;

                return map;
            }

            // Save index and patch for posterity
            slaveMap.index() = sIndex;
            slaveMap.patchIndex() = pI;
        }

        // Next collapse each slave face (for real this time...)
        forAll(slaveMaps, slaveI)
        {
            // Alias for convenience...
            changeMap& slaveMap = slaveMaps[slaveI];

            label sIndex = slaveMap.index();
            label pI = slaveMap.patchIndex();
            label slaveOverRide = slaveMap.type();

            // Temporarily turn off coupledModification
            unsetCoupledModification();

            // Collapse the slave
            if (localCouple)
            {
                slaveMap =
                (
                    collapseQuadFace(sIndex, slaveOverRide, false, forceOp)
                );
            }
            else
            {
                const coupleMap& cMap = recvMeshes_[pI].map();
                dynamicTopoFvMesh& sMesh = recvMeshes_[pI].subMesh();

                if (sIndex < 0)
                {
                    // Edge-based coupling

                    // Fetch the slave edge
                    edge sEdge = sMesh.edges_[mag(sIndex)];

                    // Move points to new location,
                    // and update operation into coupleMap
                    sMesh.points_[sEdge[0]] = slaveMoveNewPoint[slaveI][0];
                    sMesh.points_[sEdge[1]] = slaveMoveNewPoint[slaveI][1];

                    sMesh.oldPoints_[sEdge[0]] = slaveMoveOldPoint[slaveI][0];
                    sMesh.oldPoints_[sEdge[1]] = slaveMoveOldPoint[slaveI][1];

                    cMap.pushOperation
                    (
                        sEdge[0],
                        coupleMap::MOVE_POINT,
                        slaveMoveNewPoint[slaveI][0],
                        slaveMoveOldPoint[slaveI][0]
                    );

                    cMap.pushOperation
                    (
                        sEdge[1],
                        coupleMap::MOVE_POINT,
                        slaveMoveNewPoint[slaveI][1],
                        slaveMoveOldPoint[slaveI][1]
                    );

                    // Force operation to succeed
                    slaveMap.type() = 1;
                }
                else
                {
                    // Face-based coupling
                    slaveMap =
                    (
                        sMesh.collapseQuadFace
                        (
                            sIndex,
                            slaveOverRide,
                            false,
                            forceOp
                        )
                    );

                    // Push operation into coupleMap
                    switch (slaveMap.type())
                    {
                        case 1:
                        {
                            cMap.pushOperation
                            (
                                sIndex,
                                coupleMap::COLLAPSE_FIRST
                            );

                            break;
                        }

                        case 2:
                        {
                            cMap.pushOperation
                            (
                                sIndex,
                                coupleMap::COLLAPSE_SECOND
                            );

                            break;
                        }

                        case 3:
                        {
                            cMap.pushOperation
                            (
                                sIndex,
                                coupleMap::COLLAPSE_MIDPOINT
                            );

                            break;
                        }
                    }
                }
            }

            // Turn it back on.
            setCoupledModification();

            // The final operation has to succeed.
            if (slaveMap.type() <= 0)
            {
                FatalErrorIn
                (
                    "\n"
                    "const changeMap "
                    "dynamicTopoFvMesh::collapseQuadFace\n"
                    "(\n"
                    "    const label fIndex,\n"
                    "    label overRideCase,\n"
                    "    bool checkOnly,\n"
                    "    bool forceOp\n"
                    ")\n"
                )
                    << "Coupled topo-change for slave failed." << nl
                    << " Face: " << fIndex << ": " << fCheck << nl
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

    // Build a hull of cells and tri-faces that are connected to each edge
    FixedList<labelList, 2> hullCells;
    FixedList<labelList, 2> hullTriFaces;

    meshOps::constructPrismHull
    (
        checkEdgeIndex[1],
        faces_,
        cells_,
        owner_,
        neighbour_,
        edgeFaces_,
        hullTriFaces[0],
        hullCells[0]
    );

    meshOps::constructPrismHull
    (
        checkEdgeIndex[2],
        faces_,
        cells_,
        owner_,
        neighbour_,
        edgeFaces_,
        hullTriFaces[1],
        hullCells[1]
    );

    // Determine the neighbouring cells
    label c0 = owner_[fIndex], c1 = neighbour_[fIndex];

    // Define variables for the prism-face calculation
    FixedList<label,2> c0BdyIndex(-1), c0IntIndex(-1);
    FixedList<label,2> c1BdyIndex(-1), c1IntIndex(-1);
    FixedList<face,2> c0BdyFace(face(3)), c0IntFace(face(4));
    FixedList<face,2> c1BdyFace(face(3)), c1IntFace(face(4));

    // Find the prism-faces
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

    if (c1 != -1)
    {
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
    }

    // Determine if either edge belongs to a boundary
    FixedList<bool, 2> nBoundCurves(false), edgeBoundary(false);

    edgeBoundary[0] = (whichEdgePatch(checkEdgeIndex[1]) > -1);
    edgeBoundary[1] = (whichEdgePatch(checkEdgeIndex[2]) > -1);

    // Decide on collapseCase
    label collapseCase = -1;

    if (edgeBoundary[0] && !edgeBoundary[1])
    {
        collapseCase = 1;
    }
    else
    if (!edgeBoundary[0] && edgeBoundary[1])
    {
        collapseCase = 2;
    }
    else
    if (edgeBoundary[0] && edgeBoundary[1])
    {
        if (c1 != -1)
        {
            if (debug > 2)
            {
                WarningIn
                (
                    "\n"
                    "const changeMap "
                    "dynamicTopoFvMesh::collapseQuadFace\n"
                    "(\n"
                    "    const label fIndex,\n"
                    "    label overRideCase,\n"
                    "    bool checkOnly,\n"
                    "    bool forceOp\n"
                    ")\n"
                )   << "Collapsing an internal face that "
                    << "lies on two boundary patches. "
                    << "Face: " << fIndex << ": " << faces_[fIndex]
                    << endl;
            }

            // Bail out for now. If proximity based refinement is
            // switched on, mesh may be sliced at this point.
            return map;
        }

        // Check if either edge lies on a bounding curve.
        if (checkBoundingCurve(checkEdgeIndex[1]))
        {
            nBoundCurves[0] = true;
        }

        if (checkBoundingCurve(checkEdgeIndex[2]))
        {
            nBoundCurves[1] = true;
        }

        // Collapse towards a bounding curve
        if (nBoundCurves[0] && !nBoundCurves[1])
        {
            collapseCase = 1;
        }
        else
        if (!nBoundCurves[0] && nBoundCurves[1])
        {
            collapseCase = 2;
        }
        else
        if (!nBoundCurves[0] && !nBoundCurves[1])
        {
            // No bounding curves. Collapse to mid-point.
            collapseCase = 3;
        }
        else
        {
            // Two bounding curves? This might cause pinching.
            // Bail out for now.
            return map;
        }
    }
    else
    {
        // Looks like this is an interior face.
        // Collapse case [3] by default
        collapseCase = 3;
    }

    // Perform an override if requested.
    if (overRideCase != -1)
    {
        collapseCase = overRideCase;
    }

    // Configure the new point-positions
    FixedList<bool, 2> check(false);
    FixedList<FixedList<label, 2>, 2> checkPoints(FixedList<label, 2>(-1));
    FixedList<point, 2> newPoint(vector::zero);
    FixedList<point, 2> oldPoint(vector::zero);

    // Replacement check points
    FixedList<label,2> original(-1), replacement(-1);

    switch (collapseCase)
    {
        case 1:
        {
            // Collapse to the first node
            original[0] = cv2; original[1] = cv3;
            replacement[0] = cv0; replacement[1] = cv1;

            newPoint[0] = points_[replacement[0]];
            newPoint[1] = points_[replacement[1]];

            oldPoint[0] = oldPoints_[replacement[0]];
            oldPoint[1] = oldPoints_[replacement[1]];

            // Define check-points
            check[1] = true;
            checkPoints[1][0] = original[0];
            checkPoints[1][1] = original[1];

            break;
        }

        case 2:
        {
            // Collapse to the second node
            original[0] = cv0; original[1] = cv1;
            replacement[0] = cv2; replacement[1] = cv3;

            newPoint[0] = points_[replacement[0]];
            newPoint[1] = points_[replacement[1]];

            oldPoint[0] = oldPoints_[replacement[0]];
            oldPoint[1] = oldPoints_[replacement[1]];

            // Define check-points
            check[0] = true;
            checkPoints[0][0] = original[0];
            checkPoints[0][1] = original[1];

            break;
        }

        case 3:
        {
            // Collapse to the mid-point
            original[0] = cv0; original[1] = cv1;
            replacement[0] = cv2; replacement[1] = cv3;

            // Define new point-positions
            newPoint[0] =
            (
                0.5 *
                (
                    points_[original[0]]
                  + points_[replacement[0]]
                )
            );

            newPoint[1] =
            (
                0.5 *
                (
                    points_[original[1]]
                  + points_[replacement[1]]
                )
            );

            // Define old point-positions
            oldPoint[0] =
            (
                0.5 *
                (
                    oldPoints_[original[0]]
                  + oldPoints_[replacement[0]]
                )
            );

            oldPoint[1] =
            (
                0.5 *
                (
                    oldPoints_[original[1]]
                  + oldPoints_[replacement[1]]
                )
            );

            // Define check-points
            check[0] = true;
            checkPoints[0][0] = original[0];
            checkPoints[0][1] = original[1];

            check[1] = true;
            checkPoints[1][0] = replacement[0];
            checkPoints[1][1] = replacement[1];

            break;
        }

        default:
        {
            // Don't think this will ever happen.
            FatalErrorIn
            (
                "\n"
                "const changeMap "
                "dynamicTopoFvMesh::collapseQuadFace\n"
                "(\n"
                "    const label fIndex,\n"
                "    label overRideCase,\n"
                "    bool checkOnly,\n"
                "    bool forceOp\n"
                ")\n"
            )
                << "Edge: " << fIndex << ": " << faces_[fIndex]
                << ". Couldn't decide on collapseCase."
                << abort(FatalError);

            break;
        }
    }

    // Keep track of resulting cell quality,
    // if collapse is indeed feasible
    scalar collapseQuality(GREAT);

    // Check collapsibility of cells around edges
    // with the re-configured point
    forAll(check, indexI)
    {
        if (!check[indexI])
        {
            continue;
        }

        // Check whether the collapse is possible.
        if
        (
            checkCollapse
            (
                (*this),
                hullTriFaces[indexI],
                c0BdyIndex,
                c1BdyIndex,
                checkPoints[indexI],
                newPoint,
                oldPoint,
                collapseQuality,
                (c1 != -1),
                forceOp
            )
        )
        {
            map.type() = 0;
            return map;
        }
    }

    // Add a map entry of the replacements as 'addedPoints'
    //  - Used in coupled mapping
    map.addPoint(replacement[0]);
    map.addPoint(replacement[1]);

    // Are we only performing checks?
    if (checkOnly)
    {
        map.type() = collapseCase;

        if (debug > 2)
        {
            Pout<< "Face: " << fIndex
                << ":: " << faces_[fIndex] << nl
                << " collapseCase determined to be: " << collapseCase << nl
                << " Resulting quality: " << collapseQuality << nl
                << " edgeBoundary: " << edgeBoundary << nl
                << " nBoundCurves: " << nBoundCurves << nl
                << " collapsePoints: " << original << nl
                << " replacePoints: " << replacement << nl
                << endl;
        }

        return map;
    }

    if (debug > 1)
    {
        const labelList& fE = faceEdges_[fIndex];

        Pout<< nl << nl
            << "Face: " << fIndex << ": " << faces_[fIndex] << nl
            << "faceEdges: " << fE << " is to be collapsed. "
            << nl;

        Pout<< " On SubMesh: " << isSubMesh_ << nl;
        Pout<< " coupledModification: " << coupledModification_ << nl;

        label epIndex = whichPatch(fIndex);

        const polyBoundaryMesh& boundary = boundaryMesh();

        Pout<< " Patch: ";

        if (epIndex == -1)
        {
            Pout<< "Internal" << nl;
        }
        else
        if (epIndex < boundary.size())
        {
            Pout<< boundary[epIndex].name() << nl;
        }
        else
        {
            Pout<< " New patch: " << epIndex << endl;
        }

        if (debug > 2)
        {
            Pout<< nl
                << "~~~~~~~~~~~~~~~~~~~~~~~~~" << nl
                << "Hulls before modification" << nl
                << "~~~~~~~~~~~~~~~~~~~~~~~~~" << nl;

            if (debug > 3)
            {
                Pout<< "Cell [0] Boundary faces: " << nl
                    << " Face: "<< c0BdyIndex[0]
                    << "::" << c0BdyFace[0] << nl
                    << " Face: "<< c0BdyIndex[1]
                    << "::" << c0BdyFace[1] << nl;

                Pout<< "Cell [0] Interior faces: " << nl
                    << " Face: "<< c0IntIndex[0]
                    << "::" << c0IntFace[0] << nl
                    << " Face: "<< c0IntIndex[1]
                    << "::" << c0IntFace[1] << nl;

                if (c1 != -1)
                {
                    Pout<< "Cell [1] Boundary faces: " << nl
                        << " Face: "<< c1BdyIndex[0]
                        << "::" << c1BdyFace[0] << nl
                        << " Face: "<< c1BdyIndex[1]
                        << "::" << c1BdyFace[1] << nl;

                    Pout<< "Cell [1] Interior faces: " << nl
                        << " Face: "<< c1IntIndex[0]
                        << "::" << c1IntFace[0] << nl
                        << " Face: "<< c1IntIndex[1]
                        << "::" << c1IntFace[1] << nl;
                }
            }

            Pout<< nl << "Cells belonging to first Edge Hull: "
                << hullCells[0] << nl;

            forAll(hullCells[0],cellI)
            {
                const cell& firstCurCell = cells_[hullCells[0][cellI]];

                Pout<< "Cell: " << hullCells[0][cellI]
                    << ": " << firstCurCell << nl;

                forAll(firstCurCell,faceI)
                {
                    Pout<< " Face: " << firstCurCell[faceI]
                        << " : " << faces_[firstCurCell[faceI]]
                        << " fE: " << faceEdges_[firstCurCell[faceI]]
                        << nl;

                    if (debug > 3)
                    {
                        const labelList& fE = faceEdges_[firstCurCell[faceI]];

                        forAll(fE, edgeI)
                        {
                            Pout<< "  Edge: " << fE[edgeI]
                                << " : " << edges_[fE[edgeI]]
                                << nl;
                        }
                    }
                }
            }

            const labelList& firstEdgeFaces = edgeFaces_[checkEdgeIndex[1]];

            Pout<< nl << "First Edge Face Hull: "
                << firstEdgeFaces << nl;

            forAll(firstEdgeFaces,indexI)
            {
                Pout<< " Face: " << firstEdgeFaces[indexI]
                    << " : " << faces_[firstEdgeFaces[indexI]]
                    << " fE: " << faceEdges_[firstEdgeFaces[indexI]]
                    << nl;

                if (debug > 3)
                {
                    const labelList& fE = faceEdges_[firstEdgeFaces[indexI]];

                    forAll(fE, edgeI)
                    {
                        Pout<< "  Edge: " << fE[edgeI]
                            << " : " << edges_[fE[edgeI]]
                            << nl;
                    }
                }
            }

            Pout<< nl << "Cells belonging to second Edge Hull: "
                << hullCells[1] << nl;

            forAll(hullCells[1], cellI)
            {
                const cell& secondCurCell = cells_[hullCells[1][cellI]];

                Pout<< "Cell: " << hullCells[1][cellI]
                    << ": " << secondCurCell << nl;

                forAll(secondCurCell, faceI)
                {
                    Pout<< " Face: " << secondCurCell[faceI]
                        << " : " << faces_[secondCurCell[faceI]]
                        << " fE: " << faceEdges_[secondCurCell[faceI]]
                        << nl;

                    if (debug > 3)
                    {
                        const labelList& fE = faceEdges_[secondCurCell[faceI]];

                        forAll(fE, edgeI)
                        {
                            Pout<< "  Edge: " << fE[edgeI]
                                << " : " << edges_[fE[edgeI]]
                                << nl;
                        }
                    }
                }
            }

            const labelList& secondEdgeFaces = edgeFaces_[checkEdgeIndex[2]];

            Pout<< nl << "Second Edge Face Hull: "
                << secondEdgeFaces << nl;

            forAll(secondEdgeFaces, indexI)
            {
                Pout<< " Face: " << secondEdgeFaces[indexI]
                    << " : " << faces_[secondEdgeFaces[indexI]]
                    << " fE: " << faceEdges_[secondEdgeFaces[indexI]]
                    << nl;

                if (debug > 3)
                {
                    const labelList& fE = faceEdges_[secondEdgeFaces[indexI]];

                    forAll(fE, edgeI)
                    {
                        Pout<< "  Edge: " << fE[edgeI]
                            << " : " << edges_[fE[edgeI]]
                            << nl;
                    }
                }
            }

            // Write out VTK files prior to change
            if (debug > 3)
            {
                labelHashSet vtkCells;

                forAll(hullCells[0], cellI)
                {
                    if (!vtkCells.found(hullCells[0][cellI]))
                    {
                        vtkCells.insert(hullCells[0][cellI]);
                    }
                }

                forAll(hullCells[1], cellI)
                {
                    if (!vtkCells.found(hullCells[1][cellI]))
                    {
                        vtkCells.insert(hullCells[1][cellI]);
                    }
                }

                writeVTK
                (
                    Foam::name(fIndex)
                  + "_Collapse_0",
                    vtkCells.toc(),
                    3, false, true
                );
            }
        }
    }

    // Edges / Quad-faces to throw or keep during collapse
    FixedList<label,2> ends(-1);
    FixedList<label,2> faceToKeep(-1), faceToThrow(-1);
    FixedList<label,4> edgeToKeep(-1), edgeToThrow(-1);

    // Maintain a list of modified faces for mapping
    dynamicLabelList modifiedFaces(10);

    // Case 2 & 3 use identical connectivity,
    // but different point locations
    if (collapseCase == 2 || collapseCase == 3)
    {
        const labelList& firstEdgeFaces = edgeFaces_[checkEdgeIndex[1]];

        // Collapse to the second node...
        forAll(firstEdgeFaces,faceI)
        {
            // Replace point indices on faces.
            meshOps::replaceLabel
            (
                cv0,
                cv2,
                faces_[firstEdgeFaces[faceI]]
            );

            meshOps::replaceLabel
            (
                cv1,
                cv3,
                faces_[firstEdgeFaces[faceI]]
            );

            // Add an entry for mapping
            if (findIndex(modifiedFaces, firstEdgeFaces[faceI]) == -1)
            {
                modifiedFaces.append(firstEdgeFaces[faceI]);
            }

            // Determine the quad-face in cell[0] & cell[1]
            // that belongs to firstEdgeFaces
            if (firstEdgeFaces[faceI] == c0IntIndex[0])
            {
                faceToKeep[0]  = c0IntIndex[1];
                faceToThrow[0] = c0IntIndex[0];
            }

            if (firstEdgeFaces[faceI] == c0IntIndex[1])
            {
                faceToKeep[0]  = c0IntIndex[0];
                faceToThrow[0] = c0IntIndex[1];
            }

            if (c1 != -1)
            {
                if (firstEdgeFaces[faceI] == c1IntIndex[0])
                {
                    faceToKeep[1]  = c1IntIndex[1];
                    faceToThrow[1] = c1IntIndex[0];
                }

                if (firstEdgeFaces[faceI] == c1IntIndex[1])
                {
                    faceToKeep[1]  = c1IntIndex[0];
                    faceToThrow[1] = c1IntIndex[1];
                }
            }
        }

        // Find common edges between quad / tri faces...
        meshOps::findCommonEdge
        (
            c0BdyIndex[0],
            faceToKeep[0],
            faceEdges_,
            edgeToKeep[0]
        );

        meshOps::findCommonEdge
        (
            c0BdyIndex[1],
            faceToKeep[0],
            faceEdges_,
            edgeToKeep[1]
        );

        meshOps::findCommonEdge
        (
            c0BdyIndex[0],
            faceToThrow[0],
            faceEdges_,
            edgeToThrow[0]
        );

        meshOps::findCommonEdge
        (
            c0BdyIndex[1],
            faceToThrow[0],
            faceEdges_,
            edgeToThrow[1]
        );

        // Size down edgeFaces for the ends.
        meshOps::findCommonEdge
        (
            faceToThrow[0],
            faceToKeep[0],
            faceEdges_,
            ends[0]
        );

        meshOps::sizeDownList
        (
            faceToThrow[0],
            edgeFaces_[ends[0]]
        );

        if (c1 != -1)
        {
            meshOps::findCommonEdge
            (
                c1BdyIndex[0],
                faceToKeep[1],
                faceEdges_,
                edgeToKeep[2]
            );

            meshOps::findCommonEdge
            (
                c1BdyIndex[1],
                faceToKeep[1],
                faceEdges_,
                edgeToKeep[3]
            );

            meshOps::findCommonEdge
            (
                c1BdyIndex[0],
                faceToThrow[1],
                faceEdges_,
                edgeToThrow[2]
            );

            meshOps::findCommonEdge
            (
                c1BdyIndex[1],
                faceToThrow[1],
                faceEdges_,
                edgeToThrow[3]
            );

            // Size down edgeFaces for the ends.
            meshOps::findCommonEdge
            (
                faceToThrow[1],
                faceToKeep[1],
                faceEdges_,
                ends[1]
            );

            meshOps::sizeDownList
            (
                faceToThrow[1],
                edgeFaces_[ends[1]]
            );
        }

        // Correct edgeFaces for triangular faces...
        forAll(edgeToThrow, indexI)
        {
            if (edgeToThrow[indexI] == -1)
            {
                continue;
            }

            const labelList& eF = edgeFaces_[edgeToThrow[indexI]];

            label origTriFace = -1, retTriFace = -1;

            // Find the original triangular face index
            forAll(eF, faceI)
            {
                if (eF[faceI] == c0BdyIndex[0])
                {
                    origTriFace = c0BdyIndex[0];
                    break;
                }

                if (eF[faceI] == c0BdyIndex[1])
                {
                    origTriFace = c0BdyIndex[1];
                    break;
                }

                if (eF[faceI] == c1BdyIndex[0])
                {
                    origTriFace = c1BdyIndex[0];
                    break;
                }

                if (eF[faceI] == c1BdyIndex[1])
                {
                    origTriFace = c1BdyIndex[1];
                    break;
                }
            }

            // Find the retained triangular face index
            forAll(eF, faceI)
            {
                if
                (
                    (eF[faceI] != origTriFace) &&
                    (eF[faceI] != faceToThrow[0]) &&
                    (eF[faceI] != faceToThrow[1])
                )
                {
                    retTriFace = eF[faceI];
                    break;
                }
            }

            // Finally replace the face index
            if (retTriFace == -1)
            {
                // Couldn't find a retained face.
                // This must be a boundary edge, so size-down instead.
                meshOps::sizeDownList
                (
                    origTriFace,
                    edgeFaces_[edgeToKeep[indexI]]
                );
            }
            else
            {
                meshOps::replaceLabel
                (
                    origTriFace,
                    retTriFace,
                    edgeFaces_[edgeToKeep[indexI]]
                );

                meshOps::replaceLabel
                (
                    edgeToThrow[indexI],
                    edgeToKeep[indexI],
                    faceEdges_[retTriFace]
                );
            }
        }

        // Correct faceEdges / edgeFaces for quad-faces...
        forAll(firstEdgeFaces,faceI)
        {
            meshOps::replaceLabel
            (
                checkEdgeIndex[1],
                checkEdgeIndex[2],
                faceEdges_[firstEdgeFaces[faceI]]
            );

            // Renumber the edges on this face
            const labelList& fE = faceEdges_[firstEdgeFaces[faceI]];

            forAll(fE, edgeI)
            {
                if (edges_[fE[edgeI]].start() == cv0)
                {
                    edges_[fE[edgeI]].start() = cv2;
                }

                if (edges_[fE[edgeI]].end() == cv0)
                {
                    edges_[fE[edgeI]].end() = cv2;
                }

                if (edges_[fE[edgeI]].start() == cv1)
                {
                    edges_[fE[edgeI]].start() = cv3;
                }

                if (edges_[fE[edgeI]].end() == cv1)
                {
                    edges_[fE[edgeI]].end() = cv3;
                }
            }

            // Size-up edgeFaces for the replacement
            if
            (
                (firstEdgeFaces[faceI] != faceToThrow[0]) &&
                (firstEdgeFaces[faceI] != faceToThrow[1]) &&
                (firstEdgeFaces[faceI] != fIndex)
            )
            {
                meshOps::sizeUpList
                (
                    firstEdgeFaces[faceI],
                    edgeFaces_[checkEdgeIndex[2]]
                );
            }
        }

        // Remove the current face from the replacement edge
        meshOps::sizeDownList
        (
            fIndex,
            edgeFaces_[checkEdgeIndex[2]]
        );

        // Replace point labels on all triangular boundary faces.
        forAll(hullCells[0],cellI)
        {
            const cell& cellToCheck = cells_[hullCells[0][cellI]];

            forAll(cellToCheck,faceI)
            {
                face& faceToCheck = faces_[cellToCheck[faceI]];

                if (faceToCheck.size() == 3)
                {
                    forAll(faceToCheck,pointI)
                    {
                        if (faceToCheck[pointI] == cv0)
                        {
                            faceToCheck[pointI] = cv2;
                        }

                        if (faceToCheck[pointI] == cv1)
                        {
                            faceToCheck[pointI] = cv3;
                        }
                    }
                }
            }
        }

        // Now that we're done with all edges, remove them.
        removeEdge(checkEdgeIndex[0]);
        removeEdge(checkEdgeIndex[1]);
        removeEdge(checkEdgeIndex[3]);

        // Update map
        map.removeEdge(checkEdgeIndex[0]);
        map.removeEdge(checkEdgeIndex[1]);
        map.removeEdge(checkEdgeIndex[2]);

        forAll(edgeToThrow, indexI)
        {
            if (edgeToThrow[indexI] == -1)
            {
                continue;
            }

            removeEdge(edgeToThrow[indexI]);

            // Update map
            map.removeEdge(edgeToThrow[indexI]);
        }

        // Delete the two points...
        removePoint(cv0);
        removePoint(cv1);

        // Update map
        map.removePoint(cv0);
        map.removePoint(cv1);
    }
    else
    {
        const labelList& secondEdgeFaces = edgeFaces_[checkEdgeIndex[2]];

        // Collapse to the first node
        forAll(secondEdgeFaces,faceI)
        {
            // Replace point indices on faces.
            meshOps::replaceLabel
            (
                cv2,
                cv0,
                faces_[secondEdgeFaces[faceI]]
            );

            meshOps::replaceLabel
            (
                cv3,
                cv1,
                faces_[secondEdgeFaces[faceI]]
            );

            // Add an entry for mapping
            if (findIndex(modifiedFaces, secondEdgeFaces[faceI]) == -1)
            {
                modifiedFaces.append(secondEdgeFaces[faceI]);
            }

            // Determine the quad-face(s) in cell[0] & cell[1]
            // that belongs to secondEdgeFaces
            if (secondEdgeFaces[faceI] == c0IntIndex[0])
            {
                faceToKeep[0]  = c0IntIndex[1];
                faceToThrow[0] = c0IntIndex[0];
            }

            if (secondEdgeFaces[faceI] == c0IntIndex[1])
            {
                faceToKeep[0]  = c0IntIndex[0];
                faceToThrow[0] = c0IntIndex[1];
            }

            if (c1 != -1)
            {
                if (secondEdgeFaces[faceI] == c1IntIndex[0])
                {
                    faceToKeep[1]  = c1IntIndex[1];
                    faceToThrow[1] = c1IntIndex[0];
                }

                if (secondEdgeFaces[faceI] == c1IntIndex[1])
                {
                    faceToKeep[1]  = c1IntIndex[0];
                    faceToThrow[1] = c1IntIndex[1];
                }
            }
        }

        // Find common edges between quad / tri faces...
        meshOps::findCommonEdge
        (
            c0BdyIndex[0],
            faceToKeep[0],
            faceEdges_,
            edgeToKeep[0]
        );

        meshOps::findCommonEdge
        (
            c0BdyIndex[1],
            faceToKeep[0],
            faceEdges_,
            edgeToKeep[1]
        );

        meshOps::findCommonEdge
        (
            c0BdyIndex[0],
            faceToThrow[0],
            faceEdges_,
            edgeToThrow[0]
        );

        meshOps::findCommonEdge
        (
            c0BdyIndex[1],
            faceToThrow[0],
            faceEdges_,
            edgeToThrow[1]
        );

        // Size down edgeFaces for the ends.
        meshOps::findCommonEdge
        (
            faceToThrow[0],
            faceToKeep[0],
            faceEdges_,
            ends[0]
        );

        meshOps::sizeDownList
        (
            faceToThrow[0],
            edgeFaces_[ends[0]]
        );

        if (c1 != -1)
        {
            meshOps::findCommonEdge
            (
                c1BdyIndex[0],
                faceToKeep[1],
                faceEdges_,
                edgeToKeep[2]
            );

            meshOps::findCommonEdge
            (
                c1BdyIndex[1],
                faceToKeep[1],
                faceEdges_,
                edgeToKeep[3]
            );

            meshOps::findCommonEdge
            (
                c1BdyIndex[0],
                faceToThrow[1],
                faceEdges_,
                edgeToThrow[2]
            );

            meshOps::findCommonEdge
            (
                c1BdyIndex[1],
                faceToThrow[1],
                faceEdges_,
                edgeToThrow[3]
            );

            // Size down edgeFaces for the ends.
            meshOps::findCommonEdge
            (
                faceToThrow[1],
                faceToKeep[1],
                faceEdges_,
                ends[1]
            );

            meshOps::sizeDownList
            (
                faceToThrow[1],
                edgeFaces_[ends[1]]
            );
        }

        // Correct edgeFaces for triangular faces...
        forAll(edgeToThrow, indexI)
        {
            if (edgeToThrow[indexI] == -1)
            {
                continue;
            }

            const labelList& eF = edgeFaces_[edgeToThrow[indexI]];

            label origTriFace = -1, retTriFace = -1;

            // Find the original triangular face index
            forAll(eF, faceI)
            {
                if (eF[faceI] == c0BdyIndex[0])
                {
                    origTriFace = c0BdyIndex[0];
                    break;
                }

                if (eF[faceI] == c0BdyIndex[1])
                {
                    origTriFace = c0BdyIndex[1];
                    break;
                }

                if (eF[faceI] == c1BdyIndex[0])
                {
                    origTriFace = c1BdyIndex[0];
                    break;
                }

                if (eF[faceI] == c1BdyIndex[1])
                {
                    origTriFace = c1BdyIndex[1];
                    break;
                }
            }

            // Find the retained triangular face index
            forAll(eF, faceI)
            {
                if
                (
                    (eF[faceI] != origTriFace) &&
                    (eF[faceI] != faceToThrow[0]) &&
                    (eF[faceI] != faceToThrow[1])
                )
                {
                    retTriFace = eF[faceI];
                    break;
                }
            }

            // Finally replace the face/edge indices
            if (retTriFace == -1)
            {
                // Couldn't find a retained face.
                // This must be a boundary edge, so size-down instead.
                meshOps::sizeDownList
                (
                    origTriFace,
                    edgeFaces_[edgeToKeep[indexI]]
                );
            }
            else
            {
                meshOps::replaceLabel
                (
                    origTriFace,
                    retTriFace,
                    edgeFaces_[edgeToKeep[indexI]]
                );

                meshOps::replaceLabel
                (
                    edgeToThrow[indexI],
                    edgeToKeep[indexI],
                    faceEdges_[retTriFace]
                );
            }
        }

        // Correct faceEdges / edgeFaces for quad-faces...
        forAll(secondEdgeFaces,faceI)
        {
            meshOps::replaceLabel
            (
                checkEdgeIndex[2],
                checkEdgeIndex[1],
                faceEdges_[secondEdgeFaces[faceI]]
            );

            // Renumber the edges on this face
            const labelList& fE = faceEdges_[secondEdgeFaces[faceI]];

            forAll(fE, edgeI)
            {
                if (edges_[fE[edgeI]].start() == cv2)
                {
                    edges_[fE[edgeI]].start() = cv0;
                }

                if (edges_[fE[edgeI]].end() == cv2)
                {
                    edges_[fE[edgeI]].end() = cv0;
                }

                if (edges_[fE[edgeI]].start() == cv3)
                {
                    edges_[fE[edgeI]].start() = cv1;
                }

                if (edges_[fE[edgeI]].end() == cv3)
                {
                    edges_[fE[edgeI]].end() = cv1;
                }
            }

            // Size-up edgeFaces for the replacement
            if
            (
                (secondEdgeFaces[faceI] != faceToThrow[0]) &&
                (secondEdgeFaces[faceI] != faceToThrow[1]) &&
                (secondEdgeFaces[faceI] != fIndex)
            )
            {
                meshOps::sizeUpList
                (
                    secondEdgeFaces[faceI],
                    edgeFaces_[checkEdgeIndex[1]]
                );
            }
        }

        // Remove the current face from the replacement edge
        meshOps::sizeDownList
        (
            fIndex,
            edgeFaces_[checkEdgeIndex[1]]
        );

        // Replace point labels on all triangular boundary faces.
        forAll(hullCells[1], cellI)
        {
            const cell& cellToCheck = cells_[hullCells[1][cellI]];

            forAll(cellToCheck, faceI)
            {
                face& faceToCheck = faces_[cellToCheck[faceI]];

                if (faceToCheck.size() == 3)
                {
                    forAll(faceToCheck, pointI)
                    {
                        if (faceToCheck[pointI] == cv2)
                        {
                            faceToCheck[pointI] = cv0;
                        }

                        if (faceToCheck[pointI] == cv3)
                        {
                            faceToCheck[pointI] = cv1;
                        }
                    }
                }
            }
        }

        // Now that we're done with all edges, remove them.
        removeEdge(checkEdgeIndex[0]);
        removeEdge(checkEdgeIndex[2]);
        removeEdge(checkEdgeIndex[3]);

        // Update map
        map.removeEdge(checkEdgeIndex[0]);
        map.removeEdge(checkEdgeIndex[2]);
        map.removeEdge(checkEdgeIndex[3]);

        forAll(edgeToThrow, indexI)
        {
            if (edgeToThrow[indexI] == -1)
            {
                continue;
            }

            removeEdge(edgeToThrow[indexI]);

            // Update map
            map.removeEdge(edgeToThrow[indexI]);
        }

        // Delete the two points...
        removePoint(cv2);
        removePoint(cv3);

        // Update map
        map.removePoint(cv2);
        map.removePoint(cv3);
    }

    if (debug > 2)
    {
        Pout<< nl
            << "~~~~~~~~~~~~~~~~~~~~~~~~" << nl
            << "Hulls after modification" << nl
            << "~~~~~~~~~~~~~~~~~~~~~~~~" << nl;

        Pout<< nl << "Cells belonging to first Edge Hull: "
            << hullCells[0] << nl;

        forAll(hullCells[0], cellI)
        {
            const cell& firstCurCell = cells_[hullCells[0][cellI]];

            Pout<< "Cell: " << hullCells[0][cellI]
                << ": " << firstCurCell << nl;

            forAll(firstCurCell, faceI)
            {
                Pout<< " Face: " << firstCurCell[faceI]
                    << " : " << faces_[firstCurCell[faceI]]
                    << " fE: " << faceEdges_[firstCurCell[faceI]]
                    << nl;
            }
        }

        const labelList& firstEdgeFaces = edgeFaces_[checkEdgeIndex[1]];

        Pout<< nl << "First Edge Face Hull: " << firstEdgeFaces << nl;

        forAll(firstEdgeFaces, indexI)
        {
            Pout<< " Face: " << firstEdgeFaces[indexI]
                << " : " << faces_[firstEdgeFaces[indexI]]
                << " fE: " << faceEdges_[firstEdgeFaces[indexI]]
                << nl;
        }

        Pout<< nl << "Cells belonging to second Edge Hull: "
            << hullCells[1] << nl;

        forAll(hullCells[1], cellI)
        {
            const cell& secondCurCell = cells_[hullCells[1][cellI]];

            Pout<< "Cell: " << hullCells[1][cellI]
                << ": " << secondCurCell << nl;

            forAll(secondCurCell, faceI)
            {
                Pout<< " Face: " << secondCurCell[faceI]
                    << " : " << faces_[secondCurCell[faceI]]
                    << " fE: " << faceEdges_[secondCurCell[faceI]]
                    << nl;
            }
        }

        const labelList& secondEdgeFaces = edgeFaces_[checkEdgeIndex[2]];

        Pout<< nl << "Second Edge Face Hull: " << secondEdgeFaces << nl;

        forAll(secondEdgeFaces, indexI)
        {
            Pout<< " Face : " << secondEdgeFaces[indexI]
                << " : " << faces_[secondEdgeFaces[indexI]]
                << " fE: " << faceEdges_[secondEdgeFaces[indexI]]
                << nl;
        }

        Pout<< "Retained face: "
            << faceToKeep[0] << ": "
            << " owner: " << owner_[faceToKeep[0]]
            << " neighbour: " << neighbour_[faceToKeep[0]]
            << nl;

        Pout<< "Discarded face: "
            << faceToThrow[0] << ": "
            << " owner: " << owner_[faceToThrow[0]]
            << " neighbour: " << neighbour_[faceToThrow[0]]
            << nl;

        if (c1 != -1)
        {
            Pout<< "Retained face: "
                << faceToKeep[1] << ": "
                << " owner: " << owner_[faceToKeep[1]]
                << " neighbour: " << neighbour_[faceToKeep[1]]
                << nl;

            Pout<< "Discarded face: "
                << faceToThrow[1] << ": "
                << " owner: " << owner_[faceToThrow[1]]
                << " neighbour: " << neighbour_[faceToThrow[1]]
                << nl;
        }

        Pout<< endl;
    }

    // Ensure proper orientation for the two retained faces
    FixedList<label,2> cellCheck(0);

    if (owner_[faceToThrow[0]] == c0)
    {
        cellCheck[0] = neighbour_[faceToThrow[0]];

        if (owner_[faceToKeep[0]] == c0)
        {
            if
            (
                (neighbour_[faceToThrow[0]] > neighbour_[faceToKeep[0]])
             && (neighbour_[faceToKeep[0]] != -1)
            )
            {
                // This face is to be flipped
                faces_[faceToKeep[0]] =
                (
                    faces_[faceToKeep[0]].reverseFace()
                );

                owner_[faceToKeep[0]] = neighbour_[faceToKeep[0]];
                neighbour_[faceToKeep[0]] = neighbour_[faceToThrow[0]];

                setFlip(faceToKeep[0]);
            }
            else
            {
                if (neighbour_[faceToThrow[0]] != -1)
                {
                    // Keep orientation intact, and update the owner
                    owner_[faceToKeep[0]] = neighbour_[faceToThrow[0]];
                }
                else
                {
                    // This face will need to be flipped and converted
                    // to a boundary face. Flip it now, so that conversion
                    // happens later.
                    faces_[faceToKeep[0]] =
                    (
                        faces_[faceToKeep[0]].reverseFace()
                    );

                    owner_[faceToKeep[0]] = neighbour_[faceToKeep[0]];
                    neighbour_[faceToKeep[0]] = -1;

                    setFlip(faceToKeep[0]);
                }
            }
        }
        else
        {
            // Keep orientation intact, and update the neighbour
            neighbour_[faceToKeep[0]] = neighbour_[faceToThrow[0]];
        }
    }
    else
    {
        cellCheck[0] = owner_[faceToThrow[0]];

        if (neighbour_[faceToKeep[0]] == c0)
        {
            if (owner_[faceToThrow[0]] < owner_[faceToKeep[0]])
            {
                // This face is to be flipped
                faces_[faceToKeep[0]] =
                (
                    faces_[faceToKeep[0]].reverseFace()
                );

                neighbour_[faceToKeep[0]] = owner_[faceToKeep[0]];
                owner_[faceToKeep[0]] = owner_[faceToThrow[0]];

                setFlip(faceToKeep[0]);
            }
            else
            {
                // Keep orientation intact, and update the neighbour
                neighbour_[faceToKeep[0]] = owner_[faceToThrow[0]];
            }
        }
        else
        {
            // Keep orientation intact, and update the owner
            owner_[faceToKeep[0]] = owner_[faceToThrow[0]];
        }
    }

    if (c1 != -1)
    {
        if (owner_[faceToThrow[1]] == c1)
        {
            cellCheck[1] = neighbour_[faceToThrow[1]];

            if (owner_[faceToKeep[1]] == c1)
            {
                if
                (
                    (neighbour_[faceToThrow[1]] > neighbour_[faceToKeep[1]])
                 && (neighbour_[faceToKeep[1]] != -1)
                )
                {
                    // This face is to be flipped
                    faces_[faceToKeep[1]] =
                    (
                        faces_[faceToKeep[1]].reverseFace()
                    );

                    owner_[faceToKeep[1]] = neighbour_[faceToKeep[1]];
                    neighbour_[faceToKeep[1]] = neighbour_[faceToThrow[1]];

                    setFlip(faceToKeep[1]);
                }
                else
                {
                    if (neighbour_[faceToThrow[1]] != -1)
                    {
                        // Keep orientation intact, and update the owner
                        owner_[faceToKeep[1]] = neighbour_[faceToThrow[1]];
                    }
                    else
                    {
                        // This face will need to be flipped and converted
                        // to a boundary face. Flip it now, so that conversion
                        // happens later.
                        faces_[faceToKeep[1]] =
                        (
                            faces_[faceToKeep[1]].reverseFace()
                        );

                        owner_[faceToKeep[1]] = neighbour_[faceToKeep[1]];
                        neighbour_[faceToKeep[1]] = -1;

                        setFlip(faceToKeep[1]);
                    }
                }
            }
            else
            {
                // Keep orientation intact, and update the neighbour
                neighbour_[faceToKeep[1]] = neighbour_[faceToThrow[1]];
            }
        }
        else
        {
            cellCheck[1] = owner_[faceToThrow[1]];

            if (neighbour_[faceToKeep[1]] == c1)
            {
                if (owner_[faceToThrow[1]] < owner_[faceToKeep[1]])
                {
                    // This face is to be flipped
                    faces_[faceToKeep[1]] =
                    (
                        faces_[faceToKeep[1]].reverseFace()
                    );

                    neighbour_[faceToKeep[1]] = owner_[faceToKeep[1]];
                    owner_[faceToKeep[1]] = owner_[faceToThrow[1]];

                    setFlip(faceToKeep[1]);
                }
                else
                {
                    // Keep orientation intact, and update the neighbour
                    neighbour_[faceToKeep[1]] = owner_[faceToThrow[1]];
                }
            }
            else
            {
                // Keep orientation intact, and update the owner
                owner_[faceToKeep[1]] = owner_[faceToThrow[1]];
            }
        }
    }

    // Remove orphaned faces
    if (owner_[faceToKeep[0]] == -1)
    {
        const face& keepFace = faces_[faceToKeep[0]];
        const labelList& rmFE = faceEdges_[faceToKeep[0]];

        labelList keepPoints(keepFace.size(), 0);

        forAll(rmFE, edgeI)
        {
            label eIndex = rmFE[edgeI];
            labelList& eFaces = edgeFaces_[eIndex];
            const edge& checkEdge = edges_[eIndex];

            if
            (
                (eFaces.size() == 1) &&
                (eFaces[0] == faceToKeep[0])
            )
            {
                // This edge has to be removed entirely.
                removeEdge(eIndex);

                // Update map
                map.removeEdge(eIndex);
            }
            else
            {
                // Mark points
                keepPoints[keepFace.which(checkEdge[0])] = 1;
                keepPoints[keepFace.which(checkEdge[1])] = 1;

                // Size-down edgeFaces
                meshOps::sizeDownList
                (
                    faceToKeep[0],
                    eFaces
                );
            }
        }

        // Check if processor patches are being
        // converted into physical ones
        label neiProc = -1;
        label kfPatch = whichPatch(faceToKeep[0]);
        label rfPatch = whichPatch(faceToThrow[0]);

        if
        (
            (getNeighbourProcessor(kfPatch) == -1) &&
            ((neiProc = getNeighbourProcessor(rfPatch)) > -1)
        )
        {
            if (isSubMesh_)
            {
                // If this is a subMesh, and faceToKeep is on
                // a physical boundary, make a 'special' entry
                // for coupled mapping purposes.
                map.addFace
                (
                    faceToThrow[0],
                    labelList(1, (-2 - kfPatch))
                );
            }
            else
            {
                // This is a wierd overhanging cell on the master
                // processor which is being removed entirely.
                // Since the processor patch face is being converted
                // to a physical one, the slave processor needs to
                // be notified of the change.
                forAll(procIndices_, pI)
                {
                    if (procIndices_[pI] == neiProc)
                    {
                        // Find slave index from the face map
                        const coupleMap& cMap = recvMeshes_[pI].map();

                        label replaceFace =
                        (
                            cMap.findSlave
                            (
                                coupleMap::FACE,
                                faceToThrow[0]
                            )
                        );

                        if (replaceFace == -1)
                        {
                            // Something is wrong here.
                            Pout<< " Could not find slave face for: " << nl
                                << "  Master face: " << faceToThrow[0]
                                << "  :: " << faces_[faceToThrow[0]] << nl
                                << "  on proc: " << procIndices_[pI] << nl
                                << "  on converting to patch: " << kfPatch
                                << abort(FatalError);
                        }

                        cMap.pushOperation
                        (
                            replaceFace,
                            coupleMap::CONVERT_PHYSICAL,
                            kfPatch
                        );

                        break;
                    }
                }
            }
        }

        // Remove unused points
        forAll(keepPoints, pointI)
        {
            if (keepPoints[pointI] == 0)
            {
                removePoint(keepFace[pointI]);
                map.removePoint(keepFace[pointI]);
            }
        }

        // Remove the face
        removeFace(faceToKeep[0]);

        // Update map
        map.removeFace(faceToKeep[0]);
    }
    else
    if
    (
        (neighbour_[faceToKeep[0]] == -1)
     && (whichPatch(faceToKeep[0]) < 0)
    )
    {
        // Obtain a copy before adding the new face,
        // since the reference might become invalid during list resizing.
        face newFace = faces_[faceToKeep[0]];
        label newOwn = owner_[faceToKeep[0]];
        labelList newFaceEdges = faceEdges_[faceToKeep[0]];

        // This face is being converted from interior to boundary. Remove
        // from the interior list and add as a boundary face to the end.
        // Edges don't have to change, since they're all on the boundary anyway.
        label newFaceIndex =
        (
            insertFace
            (
                whichPatch(faceToThrow[0]),
                newFace,
                newOwn,
                -1,
                newFaceEdges
            )
        );

        // Add an entry for mapping
        modifiedFaces.append(newFaceIndex);

        // Update map
        map.addFace(newFaceIndex, labelList(1, faceToThrow[0]));

        meshOps::replaceLabel
        (
            faceToKeep[0],
            newFaceIndex,
            cells_[newOwn]
        );

        // Correct edgeFaces with the new face label.
        forAll(newFaceEdges, edgeI)
        {
            meshOps::replaceLabel
            (
                faceToKeep[0],
                newFaceIndex,
                edgeFaces_[newFaceEdges[edgeI]]
            );
        }

        // Renumber the neighbour so that this face is removed correctly.
        neighbour_[faceToKeep[0]] = 0;
        removeFace(faceToKeep[0]);

        // Update map
        map.removeFace(faceToKeep[0]);
    }

    // Remove the unwanted faces in the cell(s) adjacent to this face,
    // and correct the cells that contain discarded faces
    const cell& cell_0 = cells_[c0];

    forAll(cell_0,faceI)
    {
        if (cell_0[faceI] != fIndex && cell_0[faceI] != faceToKeep[0])
        {
            removeFace(cell_0[faceI]);

            // Update map
            map.removeFace(cell_0[faceI]);
        }
    }

    if (cellCheck[0] != -1)
    {
        meshOps::replaceLabel
        (
            faceToThrow[0],
            faceToKeep[0],
            cells_[cellCheck[0]]
        );
    }

    // Remove the cell
    removeCell(c0);

    // Update map
    map.removeCell(c0);

    if (c1 == -1)
    {
        // Increment the surface-collapse counter
        status(SURFACE_COLLAPSES)++;
    }
    else
    {
        // Remove orphaned faces
        if (owner_[faceToKeep[1]] == -1)
        {
            const face& keepFace = faces_[faceToKeep[1]];
            const labelList& rmFE = faceEdges_[faceToKeep[1]];

            labelList keepPoints(keepFace.size(), 0);

            forAll(rmFE, edgeI)
            {
                label eIndex = rmFE[edgeI];
                labelList& eFaces = edgeFaces_[eIndex];
                const edge& checkEdge = edges_[eIndex];

                if
                (
                    (eFaces.size() == 1) &&
                    (eFaces[0] == faceToKeep[1])
                )
                {
                    // This edge has to be removed entirely.
                    removeEdge(eIndex);

                    // Update map
                    map.removeEdge(eIndex);
                }
                else
                {
                    // Mark points
                    keepPoints[keepFace.which(checkEdge[0])] = 1;
                    keepPoints[keepFace.which(checkEdge[1])] = 1;

                    // Size-down edgeFaces
                    meshOps::sizeDownList
                    (
                        faceToKeep[1],
                        eFaces
                    );
                }
            }

            // Check if processor patches are being
            // converted into physical ones
            label neiProc = -1;
            label kfPatch = whichPatch(faceToKeep[1]);
            label rfPatch = whichPatch(faceToThrow[1]);

            if
            (
                (getNeighbourProcessor(kfPatch) == -1) &&
                ((neiProc = getNeighbourProcessor(rfPatch)) > -1)
            )
            {
                if (isSubMesh_)
                {
                    // If this is a subMesh, and faceToKeep is on
                    // a physical boundary, make a 'special' entry
                    // for coupled mapping purposes.
                    map.addFace
                    (
                        faceToThrow[1],
                        labelList(1, (-2 - kfPatch))
                    );
                }
                else
                {
                    // This is a wierd overhanging cell on the master
                    // processor which is being removed entirely.
                    // Since the processor patch face is being converted
                    // to a physical one, the slave processor needs to
                    // be notified of the change.
                    forAll(procIndices_, pI)
                    {
                        if (procIndices_[pI] == neiProc)
                        {
                            // Find slave index from the face map
                            const coupleMap& cMap = recvMeshes_[pI].map();

                            label replaceFace =
                            (
                                cMap.findSlave
                                (
                                    coupleMap::FACE,
                                    faceToThrow[1]
                                )
                            );

                            if (replaceFace == -1)
                            {
                                // Something is wrong here.
                                Pout<< " Could not find slave face for: " << nl
                                    << "  Master face: " << faceToThrow[1]
                                    << "  :: " << faces_[faceToThrow[1]] << nl
                                    << "  on proc: " << procIndices_[pI] << nl
                                    << "  on converting to patch: " << kfPatch
                                    << abort(FatalError);
                            }

                            cMap.pushOperation
                            (
                                replaceFace,
                                coupleMap::CONVERT_PHYSICAL,
                                kfPatch
                            );

                            break;
                        }
                    }
                }
            }

            // Remove unused points
            forAll(keepPoints, pointI)
            {
                if (keepPoints[pointI] == 0)
                {
                    removePoint(keepFace[pointI]);
                    map.removePoint(keepFace[pointI]);
                }
            }

            // Remove the face
            removeFace(faceToKeep[1]);

            // Update map
            map.removeFace(faceToKeep[1]);
        }
        else
        if
        (
            (neighbour_[faceToKeep[1]] == -1)
         && (whichPatch(faceToKeep[1]) < 0)
        )
        {
            // Obtain a copy before adding the new face,
            // since the reference might become invalid during list resizing.
            face newFace = faces_[faceToKeep[1]];
            label newOwn = owner_[faceToKeep[1]];
            labelList newFaceEdges = faceEdges_[faceToKeep[1]];

            // This face is being converted from interior to boundary. Remove
            // from the interior list and add as a boundary face to the end.
            // Edges don't have to change, since they're on the boundary anyway.
            label newFaceIndex =
            (
                insertFace
                (
                    whichPatch(faceToThrow[1]),
                    newFace,
                    newOwn,
                    -1,
                    newFaceEdges
                )
            );

            // Add an entry for mapping
            modifiedFaces.append(newFaceIndex);

            // Update map
            map.addFace(newFaceIndex, labelList(1, faceToThrow[1]));

            meshOps::replaceLabel
            (
                faceToKeep[1],
                newFaceIndex,
                cells_[newOwn]
            );

            // Correct edgeFaces with the new face label.
            forAll(newFaceEdges, edgeI)
            {
                meshOps::replaceLabel
                (
                    faceToKeep[1],
                    newFaceIndex,
                    edgeFaces_[newFaceEdges[edgeI]]
                );
            }

            // Renumber the neighbour so that this face is removed correctly.
            neighbour_[faceToKeep[1]] = 0;
            removeFace(faceToKeep[1]);

            // Update map
            map.removeFace(faceToKeep[1]);
        }

        const cell& cell_1 = cells_[c1];

        forAll(cell_1, faceI)
        {
            if (cell_1[faceI] != fIndex && cell_1[faceI] != faceToKeep[1])
            {
                removeFace(cell_1[faceI]);

                // Update map
                map.removeFace(cell_1[faceI]);
            }
        }

        if (cellCheck[1] != -1)
        {
            meshOps::replaceLabel
            (
                faceToThrow[1],
                faceToKeep[1],
                cells_[cellCheck[1]]
            );
        }

        // Remove the cell
        removeCell(c1);

        // Update map
        map.removeCell(c1);
    }

    // Move old / new points
    oldPoints_[replacement[0]] = oldPoint[0];
    oldPoints_[replacement[1]] = oldPoint[1];

    points_[replacement[0]] = newPoint[0];
    points_[replacement[1]] = newPoint[1];

    // Finally remove the face
    removeFace(fIndex);

    // Update map
    map.removeFace(fIndex);

    // Write out VTK files after change
    if (debug > 3)
    {
        labelHashSet vtkCells;

        forAll(hullCells[0], cellI)
        {
            if (hullCells[0][cellI] == c0 || hullCells[0][cellI] == c1)
            {
                continue;
            }

            if (!vtkCells.found(hullCells[0][cellI]))
            {
                vtkCells.insert(hullCells[0][cellI]);
            }
        }

        forAll(hullCells[1], cellI)
        {
            if (hullCells[1][cellI] == c0 || hullCells[1][cellI] == c1)
            {
                continue;
            }

            if (!vtkCells.found(hullCells[1][cellI]))
            {
                vtkCells.insert(hullCells[1][cellI]);
            }
        }

        writeVTK
        (
            Foam::name(fIndex)
          + "_Collapse_1",
            vtkCells.toc(),
            3, false, true
        );
    }

    // Fill-in candidate mapping information
    labelList mC(2, -1);
    mC[0] = c0, mC[1] = c1;

    // Now that all old / new cells possess correct connectivity,
    // compute mapping information.
    forAll(hullCells, indexI)
    {
        forAll(hullCells[indexI], cellI)
        {
            label mcIndex = hullCells[indexI][cellI];

            // Skip collapsed cells
            if (mcIndex == c0 || mcIndex == c1)
            {
                continue;
            }

            // Set the mapping for this cell
            setCellMapping(mcIndex, mC);
        }
    }

    // Set face mapping information for modified faces
    forAll(modifiedFaces, faceI)
    {
        const label mfIndex = modifiedFaces[faceI];

        // Exclude deleted faces
        if (faces_[mfIndex].empty())
        {
            continue;
        }

        // Decide between default / weighted mapping
        // based on boundary information
        label fPatch = whichPatch(mfIndex);

        if (fPatch == -1)
        {
            setFaceMapping(mfIndex);
        }
        else
        {
            // Fill-in candidate mapping information
            labelList faceCandidates;

            const labelList& fEdges = faceEdges_[mfIndex];

            forAll(fEdges, edgeI)
            {
                if (whichEdgePatch(fEdges[edgeI]) == fPatch)
                {
                    // Loop through associated edgeFaces
                    const labelList& eFaces = edgeFaces_[fEdges[edgeI]];

                    forAll(eFaces, faceI)
                    {
                        if
                        (
                            (eFaces[faceI] != mfIndex) &&
                            (whichPatch(eFaces[faceI]) == fPatch)
                        )
                        {
                            faceCandidates.setSize
                            (
                                faceCandidates.size() + 1,
                                eFaces[faceI]
                            );
                        }
                    }
                }
            }

            if (faceCandidates.empty())
            {
                // Add the face itself
                faceCandidates.setSize(1, mfIndex);
            }

            // Set the mapping for this face
            setFaceMapping(mfIndex, faceCandidates);
        }
    }

    // If modification is coupled, update mapping info.
    if (coupledModification_)
    {
        // Check if the collapse points are present
        // on a processor not involved in the current
        // operation, and update if necessary
        if (procCouple && !localCouple)
        {
            forAll(procIndices_, pI)
            {
                bool involved = false;

                forAll(slaveMaps, slaveI)
                {
                    // Alias for convenience...
                    const changeMap& slaveMap = slaveMaps[slaveI];

                    if (slaveMap.patchIndex() == pI && slaveMap.index() >= 0)
                    {
                        // Involved in this operation. Break out.
                        involved = true;
                        break;
                    }
                }

                if (involved)
                {
                    continue;
                }

                // Check coupleMaps for point coupling
                const label pointEnum = coupleMap::POINT;

                const coupledMesh& recvMesh = recvMeshes_[pI];
                const coupleMap& cMap = recvMesh.map();

                // Obtain non-const references
                Map<label>& pointMap = cMap.entityMap(pointEnum);
                Map<label>& rPointMap = cMap.reverseEntityMap(pointEnum);

                label sI0 = -1, sI1 = -1;

                if (collapsingSlave)
                {
                    if ((sI0 = cMap.findMaster(pointEnum, original[0])) > -1)
                    {
                        if (rPointMap.found(replacement[0]))
                        {
                            rPointMap[replacement[0]] = sI0;
                        }
                        else
                        {
                            rPointMap.insert(replacement[0], sI0);
                        }

                        pointMap[sI0] = replacement[0];
                    }

                    if ((sI1 = cMap.findMaster(pointEnum, original[1])) > -1)
                    {
                        if (rPointMap.found(replacement[1]))
                        {
                            rPointMap[replacement[1]] = sI1;
                        }
                        else
                        {
                            rPointMap.insert(replacement[1], sI1);
                        }

                        pointMap[sI1] = replacement[1];
                    }
                }
                else
                {
                    if ((sI0 = cMap.findSlave(pointEnum, original[0])) > -1)
                    {
                        if (pointMap.found(replacement[0]))
                        {
                            pointMap[replacement[0]] = sI0;
                        }
                        else
                        {
                            pointMap.insert(replacement[0], sI0);
                        }

                        rPointMap[sI0] = replacement[0];
                    }

                    if ((sI1 = cMap.findSlave(pointEnum, original[1])) > -1)
                    {
                        if (pointMap.found(replacement[1]))
                        {
                            pointMap[replacement[1]] = sI1;
                        }
                        else
                        {
                            pointMap.insert(replacement[1], sI1);
                        }

                        rPointMap[sI1] = replacement[1];
                    }
                }

                if (sI0 > -1 && sI1 > -1 && debug > 2)
                {
                    Pout<< " Found " << original[0] << " and " << original[1]
                        << " on proc: " << procIndices_[pI]
                        << endl;
                }
            }
        }

        forAll(slaveMaps, slaveI)
        {
            // Alias for convenience...
            const changeMap& slaveMap = slaveMaps[slaveI];

            // Skip updates for edge-based coupling
            if (slaveMap.index() < 0)
            {
                continue;
            }

            label pI = slaveMap.patchIndex();

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

            // Configure the slave replacement points.
            //  - collapseQuadFace stores this as 'addedPoints'
            label scP0 = slaveMap.removedPointList()[0];
            label scP1 = slaveMap.removedPointList()[1];

            label srP0 = slaveMap.addedPointList()[0].index();
            label srP1 = slaveMap.addedPointList()[1].index();

            // Alias for convenience
            const coupleMap& cMap = *cMapPtr;

            // Remove the point entries.
            const label pointEnum = coupleMap::POINT;

            // Obtain references
            Map<label>& pointMap = cMap.entityMap(pointEnum);
            Map<label>& rPointMap = cMap.reverseEntityMap(pointEnum);

            if (collapsingSlave)
            {
                if (rPointMap[replacement[0]] == scP0)
                {
                    pointMap[srP0] = replacement[0];
                    rPointMap[replacement[0]] = srP0;
                }
                else
                if (rPointMap[replacement[0]] == scP1)
                {
                    pointMap[srP1] = replacement[0];
                    rPointMap[replacement[0]] = srP1;
                }

                if (rPointMap[original[0]] == scP0)
                {
                    pointMap.erase(scP0);
                }
                else
                if (rPointMap[original[0]] == scP1)
                {
                    pointMap.erase(scP1);
                }

                rPointMap.erase(original[0]);
            }
            else
            {
                if (pointMap[replacement[0]] == scP0)
                {
                    rPointMap[srP0] = replacement[0];
                    pointMap[replacement[0]] = srP0;
                }
                else
                if (pointMap[replacement[0]] == scP1)
                {
                    rPointMap[srP1] = replacement[0];
                    pointMap[replacement[0]] = srP1;
                }

                if (pointMap[original[0]] == scP0)
                {
                    rPointMap.erase(scP0);
                }
                else
                if (pointMap[original[0]] == scP1)
                {
                    rPointMap.erase(scP1);
                }

                pointMap.erase(original[0]);
            }

            if (collapsingSlave)
            {
                if (rPointMap[replacement[1]] == scP0)
                {
                    pointMap[srP0] = replacement[1];
                    rPointMap[replacement[1]] = srP0;
                }
                else
                if (rPointMap[replacement[1]] == scP1)
                {
                    pointMap[srP1] = replacement[1];
                    rPointMap[replacement[1]] = srP1;
                }

                if (rPointMap[original[1]] == scP0)
                {
                    pointMap.erase(scP0);
                }
                else
                if (rPointMap[original[1]] == scP1)
                {
                    pointMap.erase(scP1);
                }

                rPointMap.erase(original[1]);
            }
            else
            {
                if (pointMap[replacement[1]] == scP0)
                {
                    rPointMap[srP0] = replacement[1];
                    pointMap[replacement[1]] = srP0;
                }
                else
                if (pointMap[replacement[1]] == scP1)
                {
                    rPointMap[srP1] = replacement[1];
                    pointMap[replacement[1]] = srP1;
                }

                if (pointMap[original[1]] == scP0)
                {
                    rPointMap.erase(scP0);
                }
                else
                if (pointMap[original[1]] == scP1)
                {
                    rPointMap.erase(scP1);
                }

                pointMap.erase(original[1]);
            }

            // Remove the face entries
            const label faceEnum = coupleMap::FACE;

            // Obtain references
            Map<label>& faceMap = cMap.entityMap(faceEnum);
            Map<label>& rFaceMap = cMap.reverseEntityMap(faceEnum);

            if (collapsingSlave)
            {
                faceMap.erase(faceMap[fIndex]);
                rFaceMap.erase(fIndex);
            }
            else
            {
                rFaceMap.erase(faceMap[fIndex]);
                faceMap.erase(fIndex);
            }

            // Configure a comparison face
            face cFace(4);

            // If any interior faces in the master map were
            // converted to boundaries, account for it
            const List<objectMap>& madF = map.addedFaceList();

            forAll(madF, faceI)
            {
                label mfIndex = madF[faceI].index();
                const face& mFace = faces_[mfIndex];

                // Select appropriate mesh
                const dynamicTopoFvMesh* meshPtr = NULL;

                // Fetch the appropriate coupleMap
                const coupleMap* crMapPtr = NULL;

                // Fetch patch info
                label ofPatch = whichPatch(fIndex);
                label mfPatch = whichPatch(mfIndex);

                if (localCouple && !procCouple)
                {
                    // Local coupling. Use this mesh itself
                    meshPtr = this;
                    crMapPtr = &(patchCoupling_[pI].map());
                }
                else
                if (procCouple && !localCouple)
                {
                    // Occasionally, this face might talk to
                    // a processor other than the slave
                    if (ofPatch == mfPatch)
                    {
                        meshPtr = &(recvMeshes_[pI].subMesh());
                        crMapPtr = &(recvMeshes_[pI].map());
                    }
                    else
                    {
                        label neiProcNo = getNeighbourProcessor(mfPatch);

                        if (neiProcNo == -1)
                        {
                            // Not a processor patch. No mapping required.
                            continue;
                        }
                        else
                        {
                            // Find an appropriate subMesh
                            label prI = findIndex(procIndices_, neiProcNo);

                            meshPtr = &(recvMeshes_[prI].subMesh());
                            crMapPtr = &(recvMeshes_[prI].map());
                        }
                    }
                }

                bool matchedFace = false;

                // Alias for convenience
                const coupleMap& crMap = *crMapPtr;
                const dynamicTopoFvMesh& sMesh = *meshPtr;

                // Obtain references
                Map<label>& pMap = crMap.entityMap(pointEnum);
                Map<label>& fMap = crMap.entityMap(faceEnum);
                Map<label>& rfMap = crMap.reverseEntityMap(faceEnum);

                // Configure the face
                forAll(cFace, pointI)
                {
                    cFace[pointI] =
                    (
                        pMap.found(mFace[pointI]) ?
                        pMap[mFace[pointI]] : -1
                    );
                }

                // Loop through all boundary faces on the subMesh
                for
                (
                    label faceJ = sMesh.nOldInternalFaces_;
                    faceJ < sMesh.faces_.size();
                    faceJ++
                )
                {
                    const face& sFace = sMesh.faces_[faceJ];

                    if (face::compare(sFace, cFace))
                    {
                        if (debug > 2)
                        {
                            Pout<< " Found face: " << faceJ
                                << "::" << sFace
                                << " with mfIndex: " << mfIndex
                                << "::" << mFace
                                << endl;
                        }

                        // Update faceMaps
                        if (collapsingSlave)
                        {
                            if (fMap.found(faceJ))
                            {
                                fMap[faceJ] = mfIndex;
                            }
                            else
                            {
                                fMap.insert(faceJ, mfIndex);
                            }

                            if (rfMap.found(mfIndex))
                            {
                                rfMap[mfIndex] = faceJ;
                            }
                            else
                            {
                                rfMap.insert(mfIndex, faceJ);
                            }
                        }
                        else
                        {
                            if (rfMap.found(faceJ))
                            {
                                rfMap[faceJ] = mfIndex;
                            }
                            else
                            {
                                rfMap.insert(faceJ, mfIndex);
                            }

                            if (fMap.found(mfIndex))
                            {
                                fMap[mfIndex] = faceJ;
                            }
                            else
                            {
                                fMap.insert(mfIndex, faceJ);
                            }
                        }

                        matchedFace = true;

                        break;
                    }
                }

                if (!matchedFace)
                {
                    // Write out for post-processing
                    writeVTK("masterFace_" + Foam::name(mfIndex), mfIndex, 2);

                    FatalErrorIn
                    (
                        "\n"
                        "const changeMap "
                        "dynamicTopoFvMesh::collapseQuadFace\n"
                        "(\n"
                        "    const label fIndex,\n"
                        "    label overRideCase,\n"
                        "    bool checkOnly,\n"
                        "    bool forceOp\n"
                        ")\n"
                    )
                        << " Master face: " << mfIndex
                        << ": " << mFace << " could not be matched." << nl
                        << " cFace: " << cFace << nl
                        << abort(FatalError);
                }
            }

            // If any interior faces in the slave map were
            // converted to boundaries, account for it
            const List<objectMap>& sadF = slaveMap.addedFaceList();

            forAll(sadF, faceI)
            {
                label sIndex = slaveMap.index();
                label sfIndex = sadF[faceI].index();

                // Select appropriate mesh
                const dynamicTopoFvMesh* meshPtr = NULL;

                if (localCouple && !procCouple)
                {
                    // Local coupling. Use this mesh itself
                    meshPtr = this;
                }
                else
                if (procCouple && !localCouple)
                {
                    meshPtr = &(recvMeshes_[pI].subMesh());
                }

                // Alias for convenience
                const dynamicTopoFvMesh& sMesh = *meshPtr;

                label ofPatch = sMesh.whichPatch(sIndex);
                label sfPatch = sMesh.whichPatch(sfIndex);

                // Skip dissimilar patches
                if (ofPatch != sfPatch && localCouple)
                {
                    continue;
                }

                const face& sFace = sMesh.faces_[sfIndex];

                if (sFace.empty())
                {
                    // Slave face was removed. Update map.
                    Map<label>::iterator sit = rFaceMap.find(sfIndex);

                    label mfIndex = -1;

                    if (sit != rFaceMap.end())
                    {
                        mfIndex = sit();
                        faceMap.erase(sit());
                        rFaceMap.erase(sit);
                    }

                    // Check if this is a special entry
                    label mo =
                    (
                        sadF[faceI].masterObjects().size() ?
                        sadF[faceI].masterObjects()[0] : 0
                    );

                    // Check if a patch conversion is necessary
                    label newPatch = -1;

                    // Back out the physical patch ID
                    if (mo < 0 && mfIndex > -1)
                    {
                        newPatch = Foam::mag(mo + 2);
                    }
                    else
                    {
                        continue;
                    }

                    // Obtain a copy before adding the new face,
                    // since the reference might become
                    // invalid during list resizing.
                    face newFace = faces_[mfIndex];
                    label newOwn = owner_[mfIndex];
                    labelList newFaceEdges = faceEdges_[mfIndex];

                    label newFaceIndex =
                    (
                        insertFace
                        (
                            newPatch,
                            newFace,
                            newOwn,
                            -1,
                            newFaceEdges
                        )
                    );

                    meshOps::replaceLabel
                    (
                        mfIndex,
                        newFaceIndex,
                        cells_[newOwn]
                    );

                    // Correct edgeFaces with the new face label.
                    forAll(newFaceEdges, edgeI)
                    {
                        meshOps::replaceLabel
                        (
                            mfIndex,
                            newFaceIndex,
                            edgeFaces_[newFaceEdges[edgeI]]
                        );
                    }

                    // Finally remove the face
                    removeFace(mfIndex);

                    // Update map
                    map.removeFace(mfIndex);
                    map.addFace(newFaceIndex, labelList(1, mfIndex));

                    continue;
                }
                else
                {
                    forAll(cFace, pointI)
                    {
                        cFace[pointI] =
                        (
                            rPointMap.found(sFace[pointI]) ?
                            rPointMap[sFace[pointI]] : -1
                        );
                    }
                }

                bool matchedFace = false;

                // Loop through all boundary faces on this mesh
                for
                (
                    label faceJ = nOldInternalFaces_;
                    faceJ < faces_.size();
                    faceJ++
                )
                {
                    const face& mFace = faces_[faceJ];

                    if (face::compare(mFace, cFace))
                    {
                        if (debug > 2)
                        {
                            Pout<< " Found face: " << faceJ
                                << "::" << mFace
                                << " with sfIndex: " << sfIndex
                                << "::" << sFace
                                << endl;
                        }

                        // Update faceMaps
                        if (collapsingSlave)
                        {
                            if (rFaceMap.found(faceJ))
                            {
                                rFaceMap[faceJ] = sfIndex;
                            }
                            else
                            {
                                rFaceMap.insert(faceJ, sfIndex);
                            }

                            if (faceMap.found(sfIndex))
                            {
                                faceMap[sfIndex] = faceJ;
                            }
                            else
                            {
                                faceMap.insert(sfIndex, faceJ);
                            }
                        }
                        else
                        {
                            if (faceMap.found(faceJ))
                            {
                                faceMap[faceJ] = sfIndex;
                            }
                            else
                            {
                                faceMap.insert(faceJ, sfIndex);
                            }

                            if (rFaceMap.found(sfIndex))
                            {
                                rFaceMap[sfIndex] = faceJ;
                            }
                            else
                            {
                                rFaceMap.insert(sfIndex, faceJ);
                            }
                        }

                        matchedFace = true;

                        break;
                    }
                }

                if (!matchedFace)
                {
                    FatalErrorIn
                    (
                        "\n"
                        "const changeMap "
                        "dynamicTopoFvMesh::collapseQuadFace\n"
                        "(\n"
                        "    const label fIndex,\n"
                        "    label overRideCase,\n"
                        "    bool checkOnly,\n"
                        "    bool forceOp\n"
                        ")\n"
                    )
                        << " Slave face: " << sfIndex
                        << ": " << sFace << " could not be matched." << nl
                        << " cFace: " << cFace << nl
                        << abort(FatalError);
                }
            }
        }
    }

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    status(TOTAL_COLLAPSES)++;

    // Increment the number of modifications
    status(TOTAL_MODIFICATIONS)++;

    // Return a succesful collapse
    map.type() = collapseCase;

    return map;
}


// Method for the collapse of an edge in 3D
// - Returns a changeMap with a type specifying:
//    -3: Collapse failed since edge was on a noRefinement patch.
//    -1: Collapse failed since max number of topo-changes was reached.
//     0: Collapse could not be performed.
//     1: Collapsed to first node.
//     2: Collapsed to second node.
//     3: Collapsed to mid-point (default).
// - overRideCase is used to force a certain collapse configuration.
//    -1: Use this value to let collapseEdge decide a case.
//     1: Force collapse to first node.
//     2: Force collapse to second node.
//     3: Force collapse to mid-point.
// - checkOnly performs a feasibility check and returns without modifications.
// - forceOp to force the collapse.
const changeMap dynamicTopoFvMesh::collapseEdge
(
    const label eIndex,
    label overRideCase,
    bool checkOnly,
    bool forceOp
)
{
    // Edge collapse performs the following operations:
    //      [1] Checks if either vertex of the edge is on a boundary
    //      [2] Checks whether cells attached to deleted vertices will be valid
    //          after the edge-collapse operation
    //      [3] Deletes all cells surrounding this edge
    //      [4] Deletes all faces surrounding this edge
    //      [5] Deletes all faces surrounding the deleted vertex attached
    //          to the cells in [3]
    //      [6] Checks the orientation of faces connected to the retained
    //          vertices
    //      [7] Remove one of the vertices of the edge
    //      Update faceEdges and edgeFaces information

    // For 2D meshes, perform face-collapse
    if (is2D())
    {
        return collapseQuadFace(eIndex, overRideCase, checkOnly);
    }

    // Figure out which thread this is...
    label tIndex = self();

    // Prepare the changeMaps
    changeMap map;
    List<changeMap> slaveMaps;
    bool collapsingSlave = false;

    if
    (
        (status(TOTAL_MODIFICATIONS) > maxModifications_)
     && (maxModifications_ > -1)
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
            "const changeMap dynamicTopoFvMesh::collapseEdge\n"
            "(\n"
            "    const label eIndex,\n"
            "    label overRideCase,\n"
            "    bool checkOnly,\n"
            "    bool forceOp\n"
            ")\n"
        )
            << " Invalid index: " << eIndex
            << abort(FatalError);
    }

    // If coupled modification is set, and this is a
    // master edge, collapse its slaves as well.
    bool localCouple = false, procCouple = false;

    if (coupledModification_)
    {
        const edge& eCheck = edges_[eIndex];

        const label edgeEnum = coupleMap::EDGE;
        const label pointEnum = coupleMap::POINT;

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

                        // Notice that we are collapsing a slave edge first.
                        collapsingSlave = true;

                        break;
                    }
                }
            }

            if (sIndex == -1)
            {
                FatalErrorIn
                (
                    "\n"
                    "const changeMap dynamicTopoFvMesh::collapseEdge\n"
                    "(\n"
                    "    const label eIndex,\n"
                    "    label overRideCase,\n"
                    "    bool checkOnly,\n"
                    "    bool forceOp\n"
                    ")\n"
                )
                    << "Coupled maps were improperly specified." << nl
                    << " Slave index not found for: " << nl
                    << " Edge: " << eIndex << ": " << eCheck << nl
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
                Pout<< nl << " >> Collapsing slave edge: " << sIndex
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
                                << "::" << eCheck
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
                else
                if
                (
                    (rcMap.findSlave(pointEnum, eCheck[0]) > -1) ||
                    (rcMap.findSlave(pointEnum, eCheck[1]) > -1)
                )
                {
                    // A point-only coupling exists.

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
                             Pout<< "Edge point on: " << eIndex
                                 << "::" << eCheck
                                 << " is handled by proc: "
                                 << procIndices_[pI]
                                 << ", so bailing out."
                                 << endl;
                         }

                         return map;
                    }

                    label p0Index = rcMap.findSlave(pointEnum, eCheck[0]);
                    label p1Index = rcMap.findSlave(pointEnum, eCheck[1]);

                    if (p0Index > -1 && p1Index == -1)
                    {
                        sIndex = p0Index;
                    }
                    else
                    if (p0Index == -1 && p1Index > -1)
                    {
                        sIndex = p1Index;
                    }

                    label curIndex = slaveMaps.size();

                    // Size up the list
                    meshOps::sizeUpList
                    (
                        changeMap(),
                        slaveMaps
                    );

                    // Save index and patch for posterity
                    //  - Negate the index to signify point coupling
                    slaveMaps[curIndex].index() = -sIndex;
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
                "const changeMap dynamicTopoFvMesh::collapseEdge\n"
                "(\n"
                "    const label eIndex,\n"
                "    label overRideCase,\n"
                "    bool checkOnly,\n"
                "    bool forceOp\n"
                ")\n"
            )
                << "Coupled maps were improperly specified." << nl
                << " localCouple: " << localCouple << nl
                << " procCouple: " << procCouple << nl
                << " Edge: " << eIndex << ": " << eCheck << nl
                << abort(FatalError);
        }

        // Temporarily turn off coupledModification
        unsetCoupledModification();

        // Test the master edge for collapse, and decide on a case
        changeMap masterMap = collapseEdge(eIndex, -1, true, forceOp);

        // Turn it back on.
        setCoupledModification();

        // Master couldn't perform collapse
        if (masterMap.type() <= 0)
        {
            return masterMap;
        }

        // For point-only coupling, define the points for checking
        pointField slaveMoveNewPoint(slaveMaps.size(), vector::zero);
        pointField slaveMoveOldPoint(slaveMaps.size(), vector::zero);

        // Now check each of the slaves for collapse feasibility
        forAll(slaveMaps, slaveI)
        {
            // Alias for convenience...
            changeMap& slaveMap = slaveMaps[slaveI];

            label slaveOverRide = -1;
            label sIndex = slaveMap.index();
            label pI = slaveMap.patchIndex();

            // If the edge is mapped onto itself, skip check
            // (occurs for cyclic edges)
            if ((sIndex == eIndex) && localCouple)
            {
                continue;
            }

            const coupleMap* cMapPtr = NULL;

            edge mEdge(eCheck), sEdge(-1, -1);

            if (localCouple)
            {
                sEdge = edges_[sIndex];

                cMapPtr = &(patchCoupling_[pI].map());
            }
            else
            if (procCouple)
            {
                const dynamicTopoFvMesh& sMesh = recvMeshes_[pI].subMesh();

                cMapPtr = &(recvMeshes_[pI].map());

                if (sIndex < 0)
                {
                    if (debug > 3)
                    {
                        Pout<< "Checking slave point: " << mag(sIndex)
                            << "::" << sMesh.points_[mag(sIndex)]
                            << " on proc: " << procIndices_[pI]
                            << " for master edge: " << eIndex
                            << " using collapseCase: " << masterMap.type()
                            << endl;
                    }
                }
                else
                {
                    sEdge = sMesh.edges_[sIndex];

                    if (debug > 3)
                    {
                        Pout<< "Checking slave edge: " << sIndex
                            << "::" << sMesh.edges_[sIndex]
                            << " on proc: " << procIndices_[pI]
                            << " for master edge: " << eIndex
                            << " using collapseCase: " << masterMap.type()
                            << endl;
                    }
                }
            }

            const Map<label>& pointMap = cMapPtr->entityMap(pointEnum);

            // Perform a topological comparison.
            switch (masterMap.type())
            {
                case 1:
                {
                    if (sIndex < 0)
                    {
                        slaveMoveNewPoint[slaveI] = points_[mEdge[0]];
                        slaveMoveOldPoint[slaveI] = oldPoints_[mEdge[0]];
                    }
                    else
                    if (pointMap[mEdge[0]] == sEdge[0])
                    {
                        slaveOverRide = 1;
                    }
                    else
                    if (pointMap[mEdge[1]] == sEdge[0])
                    {
                        slaveOverRide = 2;
                    }
                    else
                    {
                        // Write out for for post-processing
                        writeVTK("mEdge_" + Foam::name(eIndex), eIndex, 1);

                        FatalErrorIn
                        (
                            "\n"
                            "const changeMap dynamicTopoFvMesh"
                            "::collapseEdge\n"
                            "(\n"
                            "    const label eIndex,\n"
                            "    label overRideCase,\n"
                            "    bool checkOnly,\n"
                            "    bool forceOp\n"
                            ")\n"
                        )
                            << "Coupled collapse failed." << nl
                            << "Master: " << eIndex << " : " << mEdge << nl
                            << "Slave: " << sIndex << " : " << sEdge << nl
                            << abort(FatalError);
                    }

                    break;
                }

                case 2:
                {
                    if (sIndex < 0)
                    {
                        slaveMoveNewPoint[slaveI] = points_[mEdge[1]];
                        slaveMoveOldPoint[slaveI] = oldPoints_[mEdge[1]];
                    }
                    else
                    if (pointMap[mEdge[1]] == sEdge[1])
                    {
                        slaveOverRide = 2;
                    }
                    else
                    if (pointMap[mEdge[0]] == sEdge[1])
                    {
                        slaveOverRide = 1;
                    }
                    else
                    {
                        // Write out for for post-processing
                        writeVTK("mEdge_" + Foam::name(eIndex), eIndex, 1);

                        FatalErrorIn
                        (
                            "\n"
                            "const changeMap dynamicTopoFvMesh"
                            "::collapseEdge\n"
                            "(\n"
                            "    const label eIndex,\n"
                            "    label overRideCase,\n"
                            "    bool checkOnly,\n"
                            "    bool forceOp\n"
                            ")\n"
                        )
                            << "Coupled collapse failed." << nl
                            << "Master: " << eIndex << " : " << mEdge << nl
                            << "Slave: " << sIndex << " : " << sEdge << nl
                            << abort(FatalError);
                    }

                    break;
                }

                case 3:
                {
                    if (sIndex < 0)
                    {
                        slaveMoveNewPoint[slaveI] =
                        (
                            0.5 * (points_[mEdge[0]] + points_[mEdge[1]])
                        );

                        slaveMoveOldPoint[slaveI] =
                        (
                            0.5 * (oldPoints_[mEdge[0]] + oldPoints_[mEdge[1]])
                        );
                    }
                    else
                    {
                        slaveOverRide = 3;
                    }

                    break;
                }
            }

            // Temporarily turn off coupledModification
            unsetCoupledModification();

            // Test the slave edge
            if (localCouple)
            {
                slaveMap = collapseEdge(sIndex, slaveOverRide, true, forceOp);
            }
            else
            if (procCouple)
            {
                dynamicTopoFvMesh& sMesh = recvMeshes_[pI].subMesh();

                if (sIndex < 0)
                {
                    // Point-based coupling
                    scalar slaveCollapseQuality(GREAT);
                    dynamicLabelList cellsChecked(10);

                    // Check cells connected to coupled point
                    const labelList& pEdges = sMesh.pointEdges_[mag(sIndex)];

                    bool infeasible = false;

                    forAll(pEdges, edgeI)
                    {
                        const labelList& eFaces =
                        (
                            sMesh.edgeFaces_[pEdges[edgeI]]
                        );

                        // Build a list of cells to check
                        forAll(eFaces, faceI)
                        {
                            label own = sMesh.owner_[eFaces[faceI]];
                            label nei = sMesh.neighbour_[eFaces[faceI]];

                            // Check owner cell
                            if (findIndex(cellsChecked, own) == -1)
                            {
                                // Check if point movement is feasible
                                if
                                (
                                    sMesh.checkCollapse
                                    (
                                        slaveMoveNewPoint[slaveI],
                                        slaveMoveOldPoint[slaveI],
                                        mag(sIndex),
                                        own,
                                        cellsChecked,
                                        slaveCollapseQuality,
                                        forceOp
                                    )
                                )
                                {
                                    infeasible = true;
                                    break;
                                }
                            }

                            if (nei == -1)
                            {
                                continue;
                            }

                            // Check neighbour cell
                            if (findIndex(cellsChecked, nei) == -1)
                            {
                                // Check if point movement is feasible
                                if
                                (
                                    sMesh.checkCollapse
                                    (
                                        slaveMoveNewPoint[slaveI],
                                        slaveMoveOldPoint[slaveI],
                                        mag(sIndex),
                                        nei,
                                        cellsChecked,
                                        slaveCollapseQuality,
                                        forceOp
                                    )
                                )
                                {
                                    infeasible = true;
                                    break;
                                }
                            }
                        }

                        if (infeasible)
                        {
                            break;
                        }
                    }

                    if (infeasible)
                    {
                        slaveMap.type() = 0;
                    }
                    else
                    {
                        slaveMap.type() = 1;
                    }
                }
                else
                {
                    // Edge-based coupling
                    slaveMap =
                    (
                        sMesh.collapseEdge
                        (
                            sIndex,
                            slaveOverRide,
                            true,
                            forceOp
                        )
                    );
                }
            }

            // Turn it back on.
            setCoupledModification();

            if (slaveMap.type() <= 0)
            {
                // Slave couldn't perform collapse.
                map.type() = -2;

                return map;
            }

            // Save index and patch for posterity
            slaveMap.index() = sIndex;
            slaveMap.patchIndex() = pI;
        }

        // Next collapse each slave edge (for real this time...)
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

            label slaveOverRide = slaveMap.type();

            // Temporarily turn off coupledModification
            unsetCoupledModification();

            // Collapse the slave
            if (localCouple)
            {
                slaveMap = collapseEdge(sIndex, slaveOverRide, false, forceOp);
            }
            else
            {
                const coupleMap& cMap = recvMeshes_[pI].map();
                dynamicTopoFvMesh& sMesh = recvMeshes_[pI].subMesh();

                if (sIndex < 0)
                {
                    // Point based coupling

                    // Move points to new location,
                    // and update operation into coupleMap
                    sMesh.points_[mag(sIndex)] = slaveMoveNewPoint[slaveI];
                    sMesh.oldPoints_[mag(sIndex)] = slaveMoveOldPoint[slaveI];

                    cMap.pushOperation
                    (
                        mag(sIndex),
                        coupleMap::MOVE_POINT,
                        slaveMoveNewPoint[slaveI],
                        slaveMoveOldPoint[slaveI]
                    );

                    // Force operation to succeed
                    slaveMap.type() = 1;
                }
                else
                {
                    // Edge-based coupling
                    slaveMap =
                    (
                        sMesh.collapseEdge
                        (
                            sIndex,
                            slaveOverRide,
                            false,
                            forceOp
                        )
                    );

                    // Push operation into coupleMap
                    switch (slaveMap.type())
                    {
                        case 1:
                        {
                            cMap.pushOperation
                            (
                                sIndex,
                                coupleMap::COLLAPSE_FIRST
                            );

                            break;
                        }

                        case 2:
                        {
                            cMap.pushOperation
                            (
                                sIndex,
                                coupleMap::COLLAPSE_SECOND
                            );

                            break;
                        }

                        case 3:
                        {
                            cMap.pushOperation
                            (
                                sIndex,
                                coupleMap::COLLAPSE_MIDPOINT
                            );

                            break;
                        }
                    }
                }
            }

            // Turn it back on.
            setCoupledModification();

            // The final operation has to succeed.
            if (slaveMap.type() <= 0)
            {
                FatalErrorIn
                (
                    "\n"
                    "const changeMap dynamicTopoFvMesh::collapseEdge\n"
                    "(\n"
                    "    const label eIndex,\n"
                    "    label overRideCase,\n"
                    "    bool checkOnly,\n"
                    "    bool forceOp\n"
                    ")\n"
                )
                    << "Coupled topo-change for slave failed." << nl
                    << " Edge: " << eIndex << ": " << eCheck << nl
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

    // Build hullVertices for this edge
    labelList vertexHull;
    buildVertexHull(eIndex, vertexHull);

    // Hull variables
    bool found = false;
    label replaceIndex = -1, m = vertexHull.size();

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
        vertexHull,
        edgeHull,
        faceHull,
        cellHull,
        ringEntities
    );

    // Check whether points of the edge lies on a boundary
    const FixedList<bool,2> edgeBoundary = checkEdgeBoundary(eIndex);
    FixedList<label, 2> nBoundCurves(0), nProcCurves(0), checkPoints(-1);

    // Decide on collapseCase
    label collapseCase = -1;

    if (edgeBoundary[0] && !edgeBoundary[1])
    {
        collapseCase = 1;
    }
    else
    if (!edgeBoundary[0] && edgeBoundary[1])
    {
        collapseCase = 2;
    }
    else
    if (edgeBoundary[0] && edgeBoundary[1])
    {
        // If this is an interior edge with two boundary points.
        // Bail out for now. If proximity based refinement is
        // switched on, mesh may be sliced at this point.
        label edgePatch = whichEdgePatch(eIndex);

        if (edgePatch == -1)
        {
            return map;
        }

        // Check if either point lies on a bounding curve
        // Used to ensure that collapses happen towards bounding curves.
        // If the edge itself is on a bounding curve, collapse is valid.
        const edge& edgeCheck = edges_[eIndex];

        forAll(edgeCheck, pointI)
        {
            const labelList& pEdges = pointEdges_[edgeCheck[pointI]];

            forAll(pEdges, edgeI)
            {
                if
                (
                    checkBoundingCurve
                    (
                        pEdges[edgeI],
                        false,
                        &(nProcCurves[pointI])
                    )
                )
                {
                    nBoundCurves[pointI]++;
                }
            }
        }

        // Check for procCurves first
        if (!coupledModification_ && !isSubMesh_)
        {
            // Bail out for now
            if (nProcCurves[0] && nProcCurves[1])
            {
                return map;
            }

            if (nProcCurves[0] && !nProcCurves[1])
            {
                if (nBoundCurves[1])
                {
                    // Bail out for now
                    return map;
                }
                else
                {
                    collapseCase = 1;
                }
            }

            if (!nProcCurves[0] && nProcCurves[1])
            {
                if (nBoundCurves[0])
                {
                    // Bail out for now
                    return map;
                }
                else
                {
                    collapseCase = 2;
                }
            }
        }

        // If still undecided, choose based on bounding curves
        if (collapseCase == -1)
        {
            // Pick the point which is connected to more bounding curves
            if (nBoundCurves[0] > nBoundCurves[1])
            {
                collapseCase = 1;
            }
            else
            if (nBoundCurves[1] > nBoundCurves[0])
            {
                collapseCase = 2;
            }
            else
            {
                // Bounding edge: collapseEdge can collapse this edge
                collapseCase = 3;
            }
        }
    }
    else
    {
        // Looks like this is an interior edge.
        // Collapse case [3] by default
        collapseCase = 3;
    }

    // Perform an override if requested.
    if (overRideCase != -1)
    {
        collapseCase = overRideCase;
    }

    // Configure the new point-position
    point newPoint = vector::zero;
    point oldPoint = vector::zero;

    label collapsePoint = -1, replacePoint = -1;

    switch (collapseCase)
    {
        case 1:
        {
            // Collapse to the first node
            replacePoint = edges_[eIndex][0];
            collapsePoint = edges_[eIndex][1];

            newPoint = points_[replacePoint];
            oldPoint = oldPoints_[replacePoint];

            checkPoints[0] = collapsePoint;

            break;
        }

        case 2:
        {
            // Collapse to the second node
            replacePoint = edges_[eIndex][1];
            collapsePoint = edges_[eIndex][0];

            newPoint = points_[replacePoint];
            oldPoint = oldPoints_[replacePoint];

            checkPoints[0] = collapsePoint;

            break;
        }

        case 3:
        {
            // Collapse to the mid-point
            replacePoint = edges_[eIndex][1];
            collapsePoint = edges_[eIndex][0];

            newPoint =
            (
                0.5 *
                (
                    points_[replacePoint]
                  + points_[collapsePoint]
                )
            );

            oldPoint =
            (
                0.5 *
                (
                    oldPoints_[replacePoint]
                  + oldPoints_[collapsePoint]
                )
            );

            checkPoints[0] = replacePoint;
            checkPoints[1] = collapsePoint;

            break;
        }

        default:
        {
            // Don't think this will ever happen.
            FatalErrorIn
            (
                "\n"
                "const changeMap dynamicTopoFvMesh::collapseEdge\n"
                "(\n"
                "    const label eIndex,\n"
                "    label overRideCase,\n"
                "    bool checkOnly,\n"
                "    bool forceOp\n"
                ")\n"
            )
                << "Edge: " << eIndex << ": " << edges_[eIndex]
                << ". Couldn't decide on collapseCase."
                << abort(FatalError);

            break;
        }
    }

    // Loop through edges and check for feasibility of collapse
    // Also, keep track of resulting cell quality,
    // if collapse is indeed feasible
    scalar collapseQuality(GREAT);
    dynamicLabelList cellsChecked(10);

    // Add all hull cells as 'checked',
    // and therefore, feasible
    forAll(cellHull, cellI)
    {
        if (cellHull[cellI] == -1)
        {
            continue;
        }

        cellsChecked.append(cellHull[cellI]);
    }

    // Check collapsibility of cells around edges
    // with the re-configured point
    forAll(checkPoints, pointI)
    {
        if (checkPoints[pointI] == -1)
        {
            continue;
        }

        const labelList& checkPointEdges = pointEdges_[checkPoints[pointI]];

        forAll(checkPointEdges, edgeI)
        {
            const labelList& eFaces = edgeFaces_[checkPointEdges[edgeI]];

            // Build a list of cells to check
            forAll(eFaces, faceI)
            {
                label own = owner_[eFaces[faceI]];
                label nei = neighbour_[eFaces[faceI]];

                // Check owner cell
                if (findIndex(cellsChecked, own) == -1)
                {
                    // Check if a collapse is feasible
                    if
                    (
                        checkCollapse
                        (
                            newPoint,
                            oldPoint,
                            checkPoints[pointI],
                            own,
                            cellsChecked,
                            collapseQuality,
                            forceOp
                        )
                    )
                    {
                        map.type() = 0;
                        return map;
                    }
                }

                if (nei == -1)
                {
                    continue;
                }

                // Check neighbour cell
                if (findIndex(cellsChecked, nei) == -1)
                {
                    // Check if a collapse is feasible
                    if
                    (
                        checkCollapse
                        (
                            newPoint,
                            oldPoint,
                            checkPoints[pointI],
                            nei,
                            cellsChecked,
                            collapseQuality,
                            forceOp
                        )
                    )
                    {
                        map.type() = 0;
                        return map;
                    }
                }
            }
        }
    }

    // Add a map entry of the replacePoint as an 'addedPoint'
    //  - Used in coupled mapping
    map.addPoint(replacePoint);

    // Are we only performing checks?
    if (checkOnly)
    {
        // Return a succesful collapse possibility
        map.type() = collapseCase;

        // Make note of the removed point
        map.removePoint(collapsePoint);

        if (debug > 2)
        {
            Pout<< nl << "Edge: " << eIndex
                << ":: " << edges_[eIndex] << nl
                << " collapseCase determined to be: " << collapseCase << nl
                << " Resulting quality: " << collapseQuality << nl
                << " collapsePoint: " << collapsePoint << nl
                << " nBoundCurves: " << nBoundCurves << nl
                << " nProcCurves: " << nProcCurves << nl
                << endl;
        }

        return map;
    }

    // Define indices on the hull for removal / replacement
    label removeEdgeIndex = -1, replaceEdgeIndex = -1;
    label removeFaceIndex = -1, replaceFaceIndex = -1;

    if (replacePoint == edges_[eIndex][0])
    {
        replaceEdgeIndex = 0;
        replaceFaceIndex = 1;
        removeEdgeIndex = 2;
        removeFaceIndex = 3;
    }
    else
    if (replacePoint == edges_[eIndex][1])
    {
        removeEdgeIndex = 0;
        removeFaceIndex = 1;
        replaceEdgeIndex = 2;
        replaceFaceIndex = 3;
    }
    else
    {
        // Don't think this will ever happen.
        FatalErrorIn
        (
            "\n"
            "const changeMap dynamicTopoFvMesh::collapseEdge\n"
            "(\n"
            "    const label eIndex,\n"
            "    label overRideCase,\n"
            "    bool checkOnly,\n"
            "    bool forceOp\n"
            ")\n"
        )
            << "Edge: " << eIndex << ": " << edges_[eIndex]
            << ". Couldn't decide on removal / replacement indices."
            << abort(FatalError);
    }

    if (debug > 1)
    {
        Pout<< nl << nl
            << "Edge: " << eIndex << ": " << edges_[eIndex]
            << " is to be collapsed. " << nl;

        Pout<< " On SubMesh: " << isSubMesh_ << nl;
        Pout<< " coupledModification: " << coupledModification_ << nl;

        label epIndex = whichEdgePatch(eIndex);

        const polyBoundaryMesh& boundary = boundaryMesh();

        Pout<< " Patch: ";

        if (epIndex == -1)
        {
            Pout<< "Internal" << nl;
        }
        else
        if (epIndex < boundary.size())
        {
            Pout<< boundary[epIndex].name() << nl;
        }
        else
        {
            Pout<< " New patch: " << epIndex << endl;
        }

        Pout<< " nBoundCurves: " << nBoundCurves << nl
            << " nProcCurves: " << nProcCurves << nl
            << " collapseCase: " << collapseCase << nl
            << " Resulting quality: " << collapseQuality << endl;

        if (debug > 2)
        {
            Pout<< " Edges: " << edgeHull << nl
                << " Faces: " << faceHull << nl
                << " Cells: " << cellHull << nl
                << " replacePoint: " << replacePoint << nl
                << " collapsePoint: " << collapsePoint << nl
                << " checkPoints: " << checkPoints << nl;

            Pout<< " ringEntities (removed faces): " << nl;

            forAll(ringEntities[removeFaceIndex], faceI)
            {
                label fIndex = ringEntities[removeFaceIndex][faceI];

                if (fIndex == -1)
                {
                    continue;
                }

                Pout<< fIndex << ": " << faces_[fIndex] << nl;
            }

            Pout<< " ringEntities (removed edges): " << nl;

            forAll(ringEntities[removeEdgeIndex], edgeI)
            {
                label ieIndex = ringEntities[removeEdgeIndex][edgeI];

                if (ieIndex == -1)
                {
                    continue;
                }

                Pout<< ieIndex << ": " << edges_[ieIndex] << nl;
            }

            Pout<< " ringEntities (replacement faces): " << nl;

            forAll(ringEntities[replaceFaceIndex], faceI)
            {
                label fIndex = ringEntities[replaceFaceIndex][faceI];

                if (fIndex == -1)
                {
                    continue;
                }

                Pout<< fIndex << ": " << faces_[fIndex] << nl;
            }

            Pout<< " ringEntities (replacement edges): " << nl;

            forAll(ringEntities[replaceEdgeIndex], edgeI)
            {
                label ieIndex = ringEntities[replaceEdgeIndex][edgeI];

                if (ieIndex == -1)
                {
                    continue;
                }

                Pout<< ieIndex << ": " << edges_[ieIndex] << nl;
            }

            labelList& collapsePointEdges = pointEdges_[collapsePoint];

            Pout<< " pointEdges (collapsePoint): ";

            forAll(collapsePointEdges, edgeI)
            {
                Pout<< collapsePointEdges[edgeI] << " ";
            }

            Pout<< nl;

            // Write out VTK files prior to change
            //  - Using old-points for convenience in post-processing
            if (debug > 3)
            {
                writeVTK
                (
                    Foam::name(eIndex)
                  + '(' + Foam::name(collapsePoint)
                  + ',' + Foam::name(replacePoint) + ')'
                  + "_Collapse_0",
                    cellsChecked,
                    3, false, true
                );
            }
        }
    }

    if (whichEdgePatch(eIndex) > -1)
    {
        // Update number of surface collapses, if necessary.
        status(SURFACE_COLLAPSES)++;
    }

    // Maintain a list of modified faces for mapping
    dynamicLabelList modifiedFaces(10);

    // Renumber all hull faces and edges
    forAll(faceHull, indexI)
    {
        // Loop through all faces of the edge to be removed
        // and reassign them to the replacement edge
        label edgeToRemove = ringEntities[removeEdgeIndex][indexI];
        label faceToRemove = ringEntities[removeFaceIndex][indexI];
        label cellToRemove = cellHull[indexI];
        label replaceEdge = ringEntities[replaceEdgeIndex][indexI];
        label replaceFace = ringEntities[replaceFaceIndex][indexI];

        label replacePatch = whichEdgePatch(replaceEdge);
        label removePatch = whichEdgePatch(edgeToRemove);

        // Check if a patch conversion is necessary
        if (replacePatch == -1 && removePatch > -1)
        {
            if (debug > 2)
            {
                Pout<< " Edge: " << replaceEdge
                    << " :: " << edges_[replaceEdge]
                    << " is interior, but edge: " << edgeToRemove
                    << " :: " << edges_[edgeToRemove]
                    << " is on a boundary patch."
                    << endl;
            }

            // Convert patch for edge
            edge newEdge = edges_[replaceEdge];
            labelList newEdgeFaces = edgeFaces_[replaceEdge];

            // Insert the new edge
            label newEdgeIndex =
            (
                insertEdge
                (
                    removePatch,
                    newEdge,
                    newEdgeFaces
                )
            );

            // Replace faceEdges for all
            // connected faces.
            forAll(newEdgeFaces, faceI)
            {
                meshOps::replaceLabel
                (
                    replaceEdge,
                    newEdgeIndex,
                    faceEdges_[newEdgeFaces[faceI]]
                );
            }

            // Remove the edge
            removeEdge(replaceEdge);

            // Update map
            map.removeEdge(replaceEdge);
            map.addEdge(newEdgeIndex, labelList(1, replaceEdge));

            // Replace index and patch
            replaceEdge = newEdgeIndex;

            // Modify ringEntities
            ringEntities[replaceEdgeIndex][indexI] = newEdgeIndex;
        }

        const labelList& rmvEdgeFaces = edgeFaces_[edgeToRemove];

        forAll(rmvEdgeFaces, faceI)
        {
            // Replace edge labels for faces
            meshOps::replaceLabel
            (
                edgeToRemove,
                replaceEdge,
                faceEdges_[rmvEdgeFaces[faceI]]
            );

            // Loop through faces associated with this edge,
            // and renumber them as well.
            const face& faceToCheck = faces_[rmvEdgeFaces[faceI]];

            if ((replaceIndex = faceToCheck.which(collapsePoint)) > -1)
            {
                // Renumber the face...
                faces_[rmvEdgeFaces[faceI]][replaceIndex] = replacePoint;

                // Add an entry for mapping
                if (findIndex(modifiedFaces, rmvEdgeFaces[faceI]) == -1)
                {
                    modifiedFaces.append(rmvEdgeFaces[faceI]);
                }
            }

            // Hull faces should be removed for the replacement edge
            if (rmvEdgeFaces[faceI] == faceHull[indexI])
            {
                meshOps::sizeDownList
                (
                    faceHull[indexI],
                    edgeFaces_[replaceEdge]
                );

                continue;
            }

            found = false;

            // Need to avoid ring faces as well.
            forAll(ringEntities[removeFaceIndex], faceII)
            {
                if
                (
                    rmvEdgeFaces[faceI]
                 == ringEntities[removeFaceIndex][faceII]
                )
                {
                    found = true;
                    break;
                }
            }

            // Size-up the replacement edge list if the face hasn't been found.
            // These faces are connected to the edge slated for
            // removal, but do not belong to the hull.
            if (!found)
            {
                meshOps::sizeUpList
                (
                    rmvEdgeFaces[faceI],
                    edgeFaces_[replaceEdge]
                );
            }
        }

        if (cellToRemove == -1)
        {
            continue;
        }

        // Size down edgeFaces for the ring edges
        meshOps::sizeDownList
        (
            faceToRemove,
            edgeFaces_[edgeHull[indexI]]
        );

        // Ensure proper orientation of retained faces
        if (owner_[faceToRemove] == cellToRemove)
        {
            if (owner_[replaceFace] == cellToRemove)
            {
                if
                (
                    (neighbour_[faceToRemove] > neighbour_[replaceFace])
                 && (neighbour_[replaceFace] != -1)
                )
                {
                    // This face is to be flipped
                    faces_[replaceFace] = faces_[replaceFace].reverseFace();
                    owner_[replaceFace] = neighbour_[replaceFace];
                    neighbour_[replaceFace] = neighbour_[faceToRemove];

                    setFlip(replaceFace);
                }
                else
                if
                (
                    (neighbour_[faceToRemove] == -1) &&
                    (neighbour_[replaceFace] != -1)
                )
                {
                    // This interior face would need to be converted
                    // to a boundary one, and flipped as well.
                    face newFace = faces_[replaceFace].reverseFace();
                    label newOwner = neighbour_[replaceFace];
                    label newNeighbour = neighbour_[faceToRemove];
                    labelList newFE = faceEdges_[replaceFace];

                    label newFaceIndex =
                    (
                        insertFace
                        (
                            whichPatch(faceToRemove),
                            newFace,
                            newOwner,
                            newNeighbour,
                            newFE
                        )
                    );

                    // Set this face aside for mapping
                    modifiedFaces.append(newFaceIndex);

                    // Update map.
                    map.addFace(newFaceIndex, labelList(1, faceToRemove));

                    // Ensure that all edges of this face are
                    // on the boundary.
                    forAll(newFE, edgeI)
                    {
                        if (whichEdgePatch(newFE[edgeI]) == -1)
                        {
                            edge newEdge = edges_[newFE[edgeI]];
                            labelList newEF = edgeFaces_[newFE[edgeI]];

                            // Need patch information for the new edge.
                            // Find the corresponding edge in ringEntities.
                            // Note that hullEdges doesn't need to be checked,
                            // since they are common to both faces.
                            label i =
                            (
                                findIndex
                                (
                                    ringEntities[replaceEdgeIndex],
                                    newFE[edgeI]
                                )
                            );

                            label repIndex =
                            (
                                whichEdgePatch
                                (
                                    ringEntities[removeEdgeIndex][i]
                                )
                            );

                            // Insert the new edge
                            label newEdgeIndex =
                            (
                                insertEdge
                                (
                                    repIndex,
                                    newEdge,
                                    newEF
                                )
                            );

                            // Replace faceEdges for all
                            // connected faces.
                            forAll(newEF, faceI)
                            {
                                meshOps::replaceLabel
                                (
                                    newFE[edgeI],
                                    newEdgeIndex,
                                    faceEdges_[newEF[faceI]]
                                );
                            }

                            // Remove the edge
                            removeEdge(newFE[edgeI]);

                            // Update map
                            map.removeEdge(newFE[edgeI]);

                            map.addEdge
                            (
                                newEdgeIndex,
                                labelList(1, newFE[edgeI])
                            );

                            // Replace faceEdges with new edge index
                            newFE[edgeI] = newEdgeIndex;

                            // Modify ringEntities
                            ringEntities[replaceEdgeIndex][i] = newEdgeIndex;
                        }
                    }

                    // Add the new faceEdges
                    faceEdges_[newFaceIndex] = newFE;

                    // Replace edgeFaces with the new face index
                    const labelList& newFEdges = faceEdges_[newFaceIndex];

                    forAll(newFEdges, edgeI)
                    {
                        meshOps::replaceLabel
                        (
                            replaceFace,
                            newFaceIndex,
                            edgeFaces_[newFEdges[edgeI]]
                        );
                    }

                    // Remove the face
                    removeFace(replaceFace);

                    // Update map
                    map.removeFace(replaceFace);

                    // Replace label for the new owner
                    meshOps::replaceLabel
                    (
                        replaceFace,
                        newFaceIndex,
                        cells_[newOwner]
                    );

                    // Modify ringEntities and replaceFace
                    replaceFace = newFaceIndex;
                    ringEntities[replaceFaceIndex][indexI] = newFaceIndex;
                }
                else
                if
                (
                    (neighbour_[faceToRemove] == -1) &&
                    (neighbour_[replaceFace] == -1)
                )
                {
                    // Wierd overhanging cell. Since replaceFace
                    // would be an orphan if this continued, remove
                    // the face entirely.
                    labelList rmFE = faceEdges_[replaceFace];

                    forAll(rmFE, edgeI)
                    {
                        if
                        (
                            (edgeFaces_[rmFE[edgeI]].size() == 1) &&
                            (edgeFaces_[rmFE[edgeI]][0] == replaceFace)
                        )
                        {
                            // This edge has to be removed entirely.
                            removeEdge(rmFE[edgeI]);

                            // Update map
                            map.removeEdge(rmFE[edgeI]);

                            label i =
                            (
                                findIndex
                                (
                                    ringEntities[replaceEdgeIndex],
                                    rmFE[edgeI]
                                )
                            );

                            if (i > -1)
                            {
                                // Modify ringEntities
                                ringEntities[replaceEdgeIndex][i] = -1;
                            }
                        }
                        else
                        {
                            // Size-down edgeFaces
                            meshOps::sizeDownList
                            (
                                replaceFace,
                                edgeFaces_[rmFE[edgeI]]
                            );
                        }
                    }

                    // If this is a subMesh, and replaceFace is on
                    // a physical boundary, make a 'special' entry
                    // for coupled mapping purposes.
                    if (isSubMesh_)
                    {
                        label kfPatch = whichPatch(replaceFace);
                        label rfPatch = whichPatch(faceToRemove);

                        if
                        (
                            (getNeighbourProcessor(rfPatch) > -1) &&
                            (getNeighbourProcessor(kfPatch) == -1)
                        )
                        {
                            map.addFace
                            (
                                faceToRemove,
                                labelList(1, (-2 - kfPatch))
                            );
                        }
                    }

                    // Remove the face
                    removeFace(replaceFace);

                    // Update map
                    map.removeFace(replaceFace);

                    // Modify ringEntities and replaceFace
                    replaceFace = -1;
                    ringEntities[replaceFaceIndex][indexI] = -1;
                }
                else
                {
                    // Keep orientation intact, and update the owner
                    owner_[replaceFace] = neighbour_[faceToRemove];
                }
            }
            else
            if (neighbour_[faceToRemove] == -1)
            {
                // This interior face would need to be converted
                // to a boundary one, but with orientation intact.
                face newFace = faces_[replaceFace];
                label newOwner = owner_[replaceFace];
                label newNeighbour = neighbour_[faceToRemove];
                labelList newFE = faceEdges_[replaceFace];

                label newFaceIndex =
                (
                    insertFace
                    (
                        whichPatch(faceToRemove),
                        newFace,
                        newOwner,
                        newNeighbour,
                        newFE
                    )
                );

                // Set this face aside for mapping
                modifiedFaces.append(newFaceIndex);

                // Update map
                map.addFace(newFaceIndex, labelList(1, faceToRemove));

                // Ensure that all edges of this face are
                // on the boundary.
                forAll(newFE, edgeI)
                {
                    if (whichEdgePatch(newFE[edgeI]) == -1)
                    {
                        edge newEdge = edges_[newFE[edgeI]];
                        labelList newEF = edgeFaces_[newFE[edgeI]];

                        // Need patch information for the new edge.
                        // Find the corresponding edge in ringEntities.
                        // Note that hullEdges doesn't need to be checked,
                        // since they are common to both faces.
                        label i =
                        (
                            findIndex
                            (
                                ringEntities[replaceEdgeIndex],
                                newFE[edgeI]
                            )
                        );

                        label repIndex =
                        (
                            whichEdgePatch
                            (
                                ringEntities[removeEdgeIndex][i]
                            )
                        );

                        // Insert the new edge
                        label newEdgeIndex =
                        (
                            insertEdge
                            (
                                repIndex,
                                newEdge,
                                newEF
                            )
                        );

                        // Replace faceEdges for all
                        // connected faces.
                        forAll(newEF, faceI)
                        {
                            meshOps::replaceLabel
                            (
                                newFE[edgeI],
                                newEdgeIndex,
                                faceEdges_[newEF[faceI]]
                            );
                        }

                        // Remove the edge
                        removeEdge(newFE[edgeI]);

                        // Update map
                        map.removeEdge(newFE[edgeI]);

                        map.addEdge
                        (
                            newEdgeIndex,
                            labelList(1, newFE[edgeI])
                        );

                        // Replace faceEdges with new edge index
                        newFE[edgeI] = newEdgeIndex;

                        // Modify ringEntities
                        ringEntities[replaceEdgeIndex][i] = newEdgeIndex;
                    }
                }

                // Add the new faceEdges
                faceEdges_[newFaceIndex] = newFE;

                // Replace edgeFaces with the new face index
                const labelList& newFEdges = faceEdges_[newFaceIndex];

                forAll(newFEdges, edgeI)
                {
                    meshOps::replaceLabel
                    (
                        replaceFace,
                        newFaceIndex,
                        edgeFaces_[newFEdges[edgeI]]
                    );
                }

                // Remove the face
                removeFace(replaceFace);

                // Update map
                map.removeFace(replaceFace);

                // Replace label for the new owner
                meshOps::replaceLabel
                (
                    replaceFace,
                    newFaceIndex,
                    cells_[newOwner]
                );

                // Modify ringEntities and replaceFace
                replaceFace = newFaceIndex;
                ringEntities[replaceFaceIndex][indexI] = newFaceIndex;
            }
            else
            {
                // Keep orientation intact, and update the neighbour
                neighbour_[replaceFace] = neighbour_[faceToRemove];
            }

            // Update the cell
            if (neighbour_[faceToRemove] != -1)
            {
                meshOps::replaceLabel
                (
                    faceToRemove,
                    replaceFace,
                    cells_[neighbour_[faceToRemove]]
                );
            }
        }
        else
        {
            if (neighbour_[replaceFace] == cellToRemove)
            {
                if (owner_[faceToRemove] < owner_[replaceFace])
                {
                    // This face is to be flipped
                    faces_[replaceFace] = faces_[replaceFace].reverseFace();
                    neighbour_[replaceFace] = owner_[replaceFace];
                    owner_[replaceFace] = owner_[faceToRemove];

                    setFlip(replaceFace);
                }
                else
                {
                    // Keep orientation intact, and update the neighbour
                    neighbour_[replaceFace] = owner_[faceToRemove];
                }
            }
            else
            {
                // Keep orientation intact, and update the owner
                owner_[replaceFace] = owner_[faceToRemove];
            }

            // Update the cell
            meshOps::replaceLabel
            (
                faceToRemove,
                replaceFace,
                cells_[owner_[faceToRemove]]
            );
        }
    }

    // Remove all hull entities
    forAll(faceHull, indexI)
    {
        label edgeToRemove = ringEntities[removeEdgeIndex][indexI];
        label faceToRemove = ringEntities[removeFaceIndex][indexI];
        label cellToRemove = cellHull[indexI];

        if (cellToRemove != -1)
        {
            // Remove faceToRemove and associated faceEdges
            removeFace(faceToRemove);

            // Update map
            map.removeFace(faceToRemove);

            // Remove the hull cell
            removeCell(cellToRemove);

            // Update map
            map.removeCell(cellToRemove);
        }

        // Remove the hull edge and associated edgeFaces
        removeEdge(edgeToRemove);

        // Update map
        map.removeEdge(edgeToRemove);

        // Remove the hull face
        removeFace(faceHull[indexI]);

        // Update map
        map.removeFace(faceHull[indexI]);
    }

    // Loop through pointEdges for the collapsePoint,
    // and replace all occurrences with replacePoint.
    // Size-up pointEdges for the replacePoint as well.
    const labelList& pEdges = pointEdges_[collapsePoint];

    forAll(pEdges, edgeI)
    {
        // Renumber edges
        label edgeIndex = pEdges[edgeI];

        if (edgeIndex != eIndex)
        {
            if (edges_[edgeIndex][0] == collapsePoint)
            {
                edges_[edgeIndex][0] = replacePoint;

                meshOps::sizeUpList
                (
                    edgeIndex,
                    pointEdges_[replacePoint]
                );
            }
            else
            if (edges_[edgeIndex][1] == collapsePoint)
            {
                edges_[edgeIndex][1] = replacePoint;

                meshOps::sizeUpList
                (
                    edgeIndex,
                    pointEdges_[replacePoint]
                );
            }
            else
            {
                // Looks like pointEdges is inconsistent
                FatalErrorIn
                (
                    "\n"
                    "const changeMap dynamicTopoFvMesh::collapseEdge\n"
                    "(\n"
                    "    const label eIndex,\n"
                    "    label overRideCase,\n"
                    "    bool checkOnly,\n"
                    "    bool forceOp\n"
                    ")\n"
                )
                    << "pointEdges is inconsistent." << nl
                    << "Point: " << collapsePoint << nl
                    << "pointEdges: " << pEdges << nl
                    << abort(FatalError);
            }

            // Loop through faces associated with this edge,
            // and renumber them as well.
            const labelList& eFaces = edgeFaces_[edgeIndex];

            forAll(eFaces, faceI)
            {
                const face& faceToCheck = faces_[eFaces[faceI]];

                if ((replaceIndex = faceToCheck.which(collapsePoint)) > -1)
                {
                    // Renumber the face...
                    faces_[eFaces[faceI]][replaceIndex] = replacePoint;

                    // Set this face aside for mapping
                    if (findIndex(modifiedFaces, eFaces[faceI]) == -1)
                    {
                        modifiedFaces.append(eFaces[faceI]);
                    }
                }
            }
        }
    }

    // Set old / new points
    oldPoints_[replacePoint] = oldPoint;
    points_[replacePoint] = newPoint;

    // Remove the collapse point
    removePoint(collapsePoint);

    // Update map
    map.removePoint(collapsePoint);

    // Remove the edge
    removeEdge(eIndex);

    // Update map
    map.removeEdge(eIndex);

    // Check for duplicate edges connected to the replacePoint
    const labelList& rpEdges = pointEdges_[replacePoint];

    dynamicLabelList mergeFaces(10);

    forAll(rpEdges, edgeI)
    {
        const edge& eCheckI = edges_[rpEdges[edgeI]];
        const label vCheckI = eCheckI.otherVertex(replacePoint);

        for (label edgeJ = edgeI + 1; edgeJ < rpEdges.size(); edgeJ++)
        {
            const edge& eCheckJ = edges_[rpEdges[edgeJ]];
            const label vCheckJ = eCheckJ.otherVertex(replacePoint);

            if (vCheckJ == vCheckI)
            {
                labelPair efCheck(rpEdges[edgeI], rpEdges[edgeJ]);

                forAll(efCheck, edgeI)
                {
                    const labelList& eF = edgeFaces_[efCheck[edgeI]];

                    forAll(eF, faceI)
                    {
                        label patch = whichPatch(eF[faceI]);

                        if (patch == -1)
                        {
                            continue;
                        }

                        if (findIndex(mergeFaces, eF[faceI]) == -1)
                        {
                            mergeFaces.append(eF[faceI]);
                        }
                    }
                }

                if (debug > 2)
                {
                    Pout<< " Found duplicate: " << eCheckI << nl
                        << " Merge faces: " << mergeFaces << nl
                        << endl;
                }
            }
        }
    }

    // Merge faces if necessary
    if (mergeFaces.size())
    {
        mergeBoundaryFaces(mergeFaces);
    }

    // Check and remove edges with an empty edgeFaces list
    const labelList& replaceEdges = ringEntities[replaceEdgeIndex];

    forAll(replaceEdges, edgeI)
    {
        label replaceEdge = replaceEdges[edgeI];

        // If the ring edge was removed, don't bother.
        if (replaceEdge > -1)
        {
            // Account for merged edges as well
            if (edgeFaces_[replaceEdge].empty())
            {
                // Is this edge truly removed? If not, remove it.
                if (edges_[replaceEdge] != edge(-1, -1))
                {
                    if (debug > 2)
                    {
                        Pout<< " Edge: " << replaceEdge
                            << " :: " << edges_[replaceEdge]
                            << " has empty edgeFaces."
                            << endl;
                    }

                    // Remove the edge
                    removeEdge(replaceEdge);

                    // Update map
                    map.removeEdge(replaceEdge);
                }
            }
        }
    }

    // Create a list of candidates for mapping
    // using the list of checked cells
    labelList mC(cellsChecked.size());

    // For cell-mapping, exclude all hull-cells
    forAll(cellsChecked, indexI)
    {
        mC[indexI] = cellsChecked[indexI];

        if (cells_[cellsChecked[indexI]].empty())
        {
            cellsChecked[indexI] = -1;
        }
        else
        if (findIndex(cellHull, cellsChecked[indexI]) > 0)
        {
            cellsChecked[indexI] = -1;
        }
    }

    // Write out VTK files after change
    //  - Using old-points for convenience in post-processing
    if (debug > 3)
    {
        writeVTK
        (
            Foam::name(eIndex)
          + '(' + Foam::name(collapsePoint)
          + ',' + Foam::name(replacePoint) + ')'
          + "_Collapse_1",
            cellsChecked,
            3, false, true
        );
    }

    // Now that all old / new cells possess correct connectivity,
    // compute mapping information.
    forAll(cellsChecked, cellI)
    {
        if (cellsChecked[cellI] < 0)
        {
            continue;
        }

        // Set the mapping for this cell
        setCellMapping(cellsChecked[cellI], mC);
    }

    // Set face mapping information for modified faces
    forAll(modifiedFaces, faceI)
    {
        const label mfIndex = modifiedFaces[faceI];

        // Exclude deleted faces
        if (faces_[mfIndex].empty())
        {
            continue;
        }

        // Decide between default / weighted mapping
        // based on boundary information
        label fPatch = whichPatch(mfIndex);

        if (fPatch == -1)
        {
            setFaceMapping(mfIndex);
        }
        else
        {
            // Fill-in candidate mapping information
            labelList faceCandidates;

            const labelList& fEdges = faceEdges_[mfIndex];

            forAll(fEdges, edgeI)
            {
                if (whichEdgePatch(fEdges[edgeI]) == fPatch)
                {
                    // Loop through associated edgeFaces
                    const labelList& eFaces = edgeFaces_[fEdges[edgeI]];

                    forAll(eFaces, faceI)
                    {
                        if
                        (
                            (eFaces[faceI] != mfIndex) &&
                            (whichPatch(eFaces[faceI]) == fPatch)
                        )
                        {
                            faceCandidates.setSize
                            (
                                faceCandidates.size() + 1,
                                eFaces[faceI]
                            );
                        }
                    }
                }
            }

            if (faceCandidates.empty())
            {
                // Add the face itself
                faceCandidates.setSize(1, mfIndex);
            }

            // Set the mapping for this face
            setFaceMapping(mfIndex, faceCandidates);
        }
    }

    // If modification is coupled, update mapping info.
    if (coupledModification_)
    {
        // Build a list of boundary edges / faces for mapping
        dynamicLabelList checkEdges(10), checkFaces(10);

        const labelList& pEdges = pointEdges_[replacePoint];

        forAll(pEdges, edgeI)
        {
            const labelList& eFaces = edgeFaces_[pEdges[edgeI]];

            forAll(eFaces, faceI)
            {
                label fPatch = whichPatch(eFaces[faceI]);

                if (localCouple && !procCouple)
                {
                    if (!locallyCoupledEntity(eFaces[faceI], true, false, true))
                    {
                        continue;
                    }

                    // Check for cyclics
                    if (boundaryMesh()[fPatch].type() == "cyclic")
                    {
                        // Check if this is a master face
                        const coupleMap& cM = patchCoupling_[fPatch].map();
                        const Map<label>& fM = cM.entityMap(coupleMap::FACE);

                        // Only add master faces
                        // (belonging to the first half)
                        if (!fM.found(eFaces[faceI]))
                        {
                            continue;
                        }
                    }
                }
                else
                if (procCouple && !localCouple)
                {
                    if (getNeighbourProcessor(fPatch) == -1)
                    {
                        continue;
                    }
                }

                // Add face and its edges for checking
                if (findIndex(checkFaces, eFaces[faceI]) == -1)
                {
                    // Add this face
                    checkFaces.append(eFaces[faceI]);

                    const labelList& fEdges = faceEdges_[eFaces[faceI]];

                    forAll(fEdges, edgeJ)
                    {
                        if (findIndex(checkEdges, fEdges[edgeJ]) == -1)
                        {
                            checkEdges.append(fEdges[edgeJ]);
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

        if (localCouple && !procCouple)
        {
            // Alias for convenience...
            const changeMap& slaveMap = slaveMaps[0];

            const label pI = slaveMap.patchIndex();
            const coupleMap& cMap = patchCoupling_[pI].map();

            // Obtain references
            Map<label>& pointMap = cMap.entityMap(coupleMap::POINT);
            Map<label>& rPointMap = cMap.reverseEntityMap(coupleMap::POINT);

            const labelList& rpList = slaveMap.removedPointList();
            const List<objectMap>& apList = slaveMap.addedPointList();

            // Configure the slave replacement point.
            //  - collapseEdge stores this as an 'addedPoint'
            label scPoint = -1, srPoint = -1;

            if ((slaveMap.index() == eIndex) && localCouple)
            {
                // Cyclic edge at axis
                scPoint = collapsePoint;
                srPoint = replacePoint;
            }
            else
            {
                scPoint = rpList[0];
                srPoint = apList[0].index();
            }

            if (collapsingSlave)
            {
                if (rPointMap[replacePoint] == scPoint)
                {
                    pointMap[srPoint] = replacePoint;
                    rPointMap[replacePoint] = srPoint;
                }
            }
            else
            {
                if (pointMap[replacePoint] == scPoint)
                {
                    rPointMap[srPoint] = replacePoint;
                    pointMap[replacePoint] = srPoint;
                }
            }

            // Update face maps
            const label faceEnum = coupleMap::FACE;

            // Obtain references
            Map<label>& faceMap = cMap.entityMap(faceEnum);
            Map<label>& rFaceMap = cMap.reverseEntityMap(faceEnum);

            forAll(checkFaces, faceI)
            {
                label mfIndex = checkFaces[faceI];
                label mfPatch = whichPatch(mfIndex);

                const face& mF = faces_[mfIndex];

                triFace cF
                (
                    pointMap.found(mF[0]) ? pointMap[mF[0]] : -1,
                    pointMap.found(mF[1]) ? pointMap[mF[1]] : -1,
                    pointMap.found(mF[2]) ? pointMap[mF[2]] : -1
                );

                if (cF[0] == -1 || cF[1] == -1 || cF[2] == -1)
                {
                    writeVTK
                    (
                        "failedFace_"
                      + Foam::name(mfIndex),
                        mfIndex,
                        2, false, true
                    );

                    Pout<< " Failed to configure face for: "
                        << mfIndex << " :: " << faces_[mfIndex]
                        << " using comparison face: " << cF
                        << abort(FatalError);
                }

                // Fetch edges connected to first slave point
                const labelList& spEdges = pointEdges_[cF[0]];

                forAll(spEdges, edgeJ)
                {
                    label seIndex = spEdges[edgeJ];

                    if (whichEdgePatch(seIndex) == -1)
                    {
                        continue;
                    }

                    const labelList& seFaces = edgeFaces_[seIndex];

                    forAll(seFaces, faceJ)
                    {
                        label sfIndex = seFaces[faceJ];

                        if (whichPatch(sfIndex) == -1)
                        {
                            continue;
                        }

                        const face& sF = faces_[sfIndex];

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
                                    << " :: " << faces_[mfIndex]
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

                if (!matchedFaces[faceI])
                {
                    writeVTK
                    (
                        "failedFacePoints_"
                      + Foam::name(mfIndex),
                        labelList(cF), 0, false, true
                    );

                    Pout<< " Failed to match face: "
                        << mfIndex << " :: " << faces_[mfIndex]
                        << " using comparison face: " << cF
                        << abort(FatalError);
                }
            }

            // Update edge maps
            const label edgeEnum = coupleMap::EDGE;

            // Obtain references
            Map<label>& edgeMap = cMap.entityMap(edgeEnum);
            Map<label>& rEdgeMap = cMap.reverseEntityMap(edgeEnum);

            forAll(checkEdges, edgeI)
            {
                label meIndex = checkEdges[edgeI];

                const edge& mE = edges_[meIndex];

                edge cE
                (
                    pointMap.found(mE[0]) ? pointMap[mE[0]] : -1,
                    pointMap.found(mE[1]) ? pointMap[mE[1]] : -1
                );

                if (cE[0] == -1 || cE[1] == -1)
                {
                    writeVTK
                    (
                        "failedEdge_"
                      + Foam::name(meIndex),
                        meIndex,
                        1, false, true
                    );

                    Pout<< " Failed to configure edge for: "
                        << meIndex << " :: " << edges_[meIndex]
                        << " using comparison edge: " << cE
                        << abort(FatalError);
                }

                // Fetch edges connected to first slave point
                const labelList& spEdges = pointEdges_[cE[0]];

                forAll(spEdges, edgeJ)
                {
                    label seIndex = spEdges[edgeJ];

                    const edge& sE = edges_[seIndex];

                    if (sE == cE)
                    {
                        if (debug > 2)
                        {
                            Pout<< " Found edge: " << seIndex
                                << " :: " << sE
                                << " with meIndex: " << meIndex
                                << " :: " << edges_[meIndex]
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
                    writeVTK
                    (
                        "failedEdgePoints_"
                      + Foam::name(meIndex),
                        labelList(cE), 0, false, true
                    );

                    Pout<< " Failed to match edge: "
                        << meIndex << " :: " << edges_[meIndex]
                        << " using comparison edge: " << cE
                        << abort(FatalError);
                }
            }
        }
        else
        if (procCouple && !localCouple)
        {
            // Update point mapping
            forAll(procIndices_, pI)
            {
                const coupleMap& cMap = recvMeshes_[pI].map();

                // Obtain references
                Map<label>& pointMap = cMap.entityMap(coupleMap::POINT);
                Map<label>& rPointMap = cMap.reverseEntityMap(coupleMap::POINT);

                const changeMap* slaveMapPtr = NULL;
                const label pointEnum = coupleMap::POINT;

                forAll(slaveMaps, slaveI)
                {
                    const changeMap& slaveMap = slaveMaps[slaveI];

                    if (slaveMap.patchIndex() == pI)
                    {
                        if (slaveMap.index() < 0)
                        {
                            // Point-coupling
                            label sI = -1;

                            if (collapsingSlave)
                            {
                                sI = cMap.findMaster(pointEnum, collapsePoint);

                                if (sI > -1)
                                {
                                    if (rPointMap.found(replacePoint))
                                    {
                                        rPointMap[replacePoint] = sI;
                                    }
                                    else
                                    {
                                        rPointMap.insert(replacePoint, sI);
                                    }

                                    pointMap[sI] = replacePoint;
                                }
                            }
                            else
                            {
                                sI = cMap.findSlave(pointEnum, collapsePoint);

                                if (sI > -1)
                                {
                                    if (pointMap.found(replacePoint))
                                    {
                                        pointMap[replacePoint] = sI;
                                    }
                                    else
                                    {
                                        pointMap.insert(replacePoint, sI);
                                    }

                                    rPointMap[sI] = replacePoint;
                                }
                            }

                            if (sI > -1 && debug > 2)
                            {
                                Pout<< " Found point: " << collapsePoint
                                    << " on proc: " << procIndices_[pI]
                                    << endl;
                            }
                        }
                        else
                        {
                            // Edge-coupling. Fetch address for later.
                            slaveMapPtr = &slaveMap;
                            break;
                        }
                    }
                }

                if (slaveMapPtr)
                {
                    // Alias for convenience...
                    const changeMap& slaveMap = *slaveMapPtr;

                    const labelList& rpList = slaveMap.removedPointList();
                    const List<objectMap>& apList = slaveMap.addedPointList();

                    // Configure the slave replacement point.
                    //  - collapseEdge stores this as an 'addedPoint'
                    label scPoint = rpList[0];
                    label srPoint = apList[0].index();

                    if (collapsingSlave)
                    {
                        if (rPointMap[replacePoint] == scPoint)
                        {
                            pointMap[srPoint] = replacePoint;
                            rPointMap[replacePoint] = srPoint;
                        }

                        pointMap.erase(rPointMap[collapsePoint]);
                        rPointMap.erase(collapsePoint);
                    }
                    else
                    {
                        if (pointMap[replacePoint] == scPoint)
                        {
                            rPointMap[srPoint] = replacePoint;
                            pointMap[replacePoint] = srPoint;
                        }

                        rPointMap.erase(pointMap[collapsePoint]);
                        pointMap.erase(collapsePoint);
                    }

                    // If any other points were removed, update map
                    for (label pointI = 1; pointI < rpList.size(); pointI++)
                    {
                        if (collapsingSlave)
                        {
                            if (pointMap.found(rpList[pointI]))
                            {
                                rPointMap.erase(pointMap[rpList[pointI]]);
                                pointMap.erase(rpList[pointI]);
                            }
                        }
                        else
                        {
                            if (rPointMap.found(rpList[pointI]))
                            {
                                if (debug > 2)
                                {
                                    Pout<< " Found removed point: "
                                        << rpList[pointI]
                                        << " on proc: " << procIndices_[pI]
                                        << " for point on this proc: "
                                        << rPointMap[rpList[pointI]]
                                        << endl;
                                }

                                pointMap.erase(rPointMap[rpList[pointI]]);
                                rPointMap.erase(rpList[pointI]);
                            }
                        }
                    }
                }
            }

            // Update face mapping
            forAll(procIndices_, pI)
            {
                const coupleMap& cMap = recvMeshes_[pI].map();
                const dynamicTopoFvMesh& sMesh = recvMeshes_[pI].subMesh();

                // Obtain point maps
                Map<label>& pointMap = cMap.entityMap(coupleMap::POINT);

                // Update face maps
                const label faceEnum = coupleMap::FACE;

                // Obtain references
                Map<label>& faceMap = cMap.entityMap(faceEnum);
                Map<label>& rFaceMap = cMap.reverseEntityMap(faceEnum);

                const changeMap* slaveMapPtr = NULL;

                forAll(slaveMaps, slaveI)
                {
                    const changeMap& slaveMap = slaveMaps[slaveI];

                    if (slaveMap.patchIndex() == pI)
                    {
                        if (slaveMap.index() < 0)
                        {
                            // Point-coupling
                            continue;
                        }

                        // Edge-coupling. Fetch address for later.
                        slaveMapPtr = &slaveMap;
                        break;
                    }
                }

                forAll(checkFaces, faceI)
                {
                    label mfIndex = checkFaces[faceI];

                    const face& mF = faces_[mfIndex];

                    // Skip checks if a conversion occurred,
                    // and this face was removed as a result
                    if (mF.empty())
                    {
                        continue;
                    }

                    label mfPatch = whichPatch(mfIndex);
                    label neiProc = getNeighbourProcessor(mfPatch);

                    triFace cF
                    (
                        pointMap.found(mF[0]) ? pointMap[mF[0]] : -1,
                        pointMap.found(mF[1]) ? pointMap[mF[1]] : -1,
                        pointMap.found(mF[2]) ? pointMap[mF[2]] : -1
                    );

                    // Check if a patch conversion is necessary
                    label newPatch = -1;
                    bool requireConversion = false, physConvert = false;

                    // Skip mapping if all points were not found
                    if (cF[0] == -1 || cF[1] == -1 || cF[2] == -1)
                    {
                        // Is this actually a non-processor patch?
                        if (neiProc == -1)
                        {
                            continue;
                        }
                        else
                        {
                            physConvert = true;
                        }
                    }

                    // Has this face been converted to a physical boundary?
                    if (faceMap.found(mfIndex) && slaveMapPtr)
                    {
                        // Alias for convenience...
                        const changeMap& sMap = *slaveMapPtr;
                        const labelList& rfList = sMap.removedFaceList();
                        const List<objectMap>& afList = sMap.addedFaceList();

                        // Was the face removed by the slave?
                        label sfIndex = faceMap[mfIndex];

                        if (findIndex(rfList, sfIndex) > -1)
                        {
                            if (debug > 2)
                            {
                                Pout<< " Found removed face: " << sfIndex
                                    << " with mfIndex: " << mfIndex
                                    << " :: " << faces_[mfIndex]
                                    << " on proc: " << procIndices_[pI]
                                    << endl;
                            }

                            // Remove map entry
                            rFaceMap.erase(sfIndex);
                            faceMap.erase(mfIndex);
                        }

                        // Check added faces for special entries
                        forAll(afList, indexI)
                        {
                            const objectMap& af = afList[indexI];

                            label mo =
                            (
                                af.masterObjects().size() ?
                                af.masterObjects()[0] : 0
                            );

                            // Back out the physical patch ID
                            if ((af.index() == sfIndex) && (mo < 0))
                            {
                                newPatch = Foam::mag(mo + 2);
                                requireConversion = true;
                                break;
                            }
                        }
                    }

                    // Was an appropriate physical patch found?
                    if (physConvert && !requireConversion)
                    {
                        continue;
                    }

                    // Are we talking to a different processor?
                    if (neiProc != procIndices_[pI])
                    {
                        // This face needs to be converted
                        const polyBoundaryMesh& boundary = boundaryMesh();

                        forAll(boundary, pi)
                        {
                            if (getNeighbourProcessor(pi) == procIndices_[pI])
                            {
                                newPatch = pi;
                                requireConversion = true;
                                break;
                            }
                        }
                    }

                    if (requireConversion)
                    {
                        // Obtain a copy before adding the new face,
                        // since the reference might become
                        // invalid during list resizing.
                        face newFace = faces_[mfIndex];
                        label newOwn = owner_[mfIndex];
                        labelList newFaceEdges = faceEdges_[mfIndex];

                        label newFaceIndex =
                        (
                            insertFace
                            (
                                newPatch,
                                newFace,
                                newOwn,
                                -1,
                                newFaceEdges
                            )
                        );

                        meshOps::replaceLabel
                        (
                            mfIndex,
                            newFaceIndex,
                            cells_[newOwn]
                        );

                        // Correct edgeFaces with the new face label.
                        forAll(newFaceEdges, edgeI)
                        {
                            meshOps::replaceLabel
                            (
                                mfIndex,
                                newFaceIndex,
                                edgeFaces_[newFaceEdges[edgeI]]
                            );
                        }

                        // Finally remove the face
                        removeFace(mfIndex);

                        // Update map
                        map.removeFace(mfIndex);
                        map.addFace(newFaceIndex, labelList(1, mfIndex));

                        // Replace index and patch
                        mfIndex = newFaceIndex;
                        mfPatch = newPatch;

                        // If conversion was to a physical patch,
                        // skip the remaining face mapping steps
                        if (getNeighbourProcessor(newPatch) == -1)
                        {
                            continue;
                        }

                        // Fail for now
                        Pout<< " Conversion required... "
                            << " Face: " << newFace << " : "
                            << newFace.centre(points_)
                            << abort(FatalError);

                        // Push a patch-conversion operation
                        cMap.pushOperation
                        (
                            newFaceIndex,
                            coupleMap::CONVERT_PATCH,
                            newFace.centre(points_),
                            newFace.centre(oldPoints_)
                        );
                    }

                    // Fetch edges connected to first slave point
                    const labelList& spEdges = sMesh.pointEdges_[cF[0]];

                    forAll(spEdges, edgeJ)
                    {
                        label seIndex = spEdges[edgeJ];

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
                                        << " :: " << faces_[mfIndex]
                                        << " on proc: " << procIndices_[pI]
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

                        Pout<< " Failed to match face: "
                            << mfIndex << " :: " << faces_[mfIndex]
                            << " using comparison face: " << cF
                            << " on proc: " << procIndices_[pI]
                            << endl;
                    }
                }
            }

            // Update edge mapping
            forAll(procIndices_, pI)
            {
                const coupleMap& cMap = recvMeshes_[pI].map();
                const dynamicTopoFvMesh& sMesh = recvMeshes_[pI].subMesh();

                // Obtain point maps
                Map<label>& pointMap = cMap.entityMap(coupleMap::POINT);

                // Update edge maps
                const label edgeEnum = coupleMap::EDGE;

                // Obtain references
                Map<label>& edgeMap = cMap.entityMap(edgeEnum);
                Map<label>& rEdgeMap = cMap.reverseEntityMap(edgeEnum);

                const changeMap* slaveMapPtr = NULL;

                forAll(slaveMaps, slaveI)
                {
                    const changeMap& slaveMap = slaveMaps[slaveI];

                    if (slaveMap.patchIndex() == pI)
                    {
                        if (slaveMap.index() < 0)
                        {
                            // Point-coupling
                            continue;
                        }

                        // Edge-coupling. Fetch address for later.
                        slaveMapPtr = &slaveMap;
                        break;
                    }
                }

                forAll(checkEdges, edgeI)
                {
                    label meIndex = checkEdges[edgeI];

                    const edge& mE = edges_[meIndex];

                    // Skip checks if a conversion occurred,
                    // and this edge was removed as a result
                    if (edgeFaces_[meIndex].empty())
                    {
                        continue;
                    }

                    label mePatch = whichEdgePatch(meIndex);
                    label neiProc = getNeighbourProcessor(mePatch);

                    edge cE
                    (
                        pointMap.found(mE[0]) ? pointMap[mE[0]] : -1,
                        pointMap.found(mE[1]) ? pointMap[mE[1]] : -1
                    );

                    // Check if a patch conversion is necessary
                    label newPatch = -1;
                    bool requireConversion = true, physConvert = false;

                    // Skip mapping if all points were not found
                    if (cE[0] == -1 || cE[1] == -1)
                    {
                        // Is this actually a non-processor patch?
                        if (neiProc == -1)
                        {
                            continue;
                        }
                        else
                        {
                            physConvert = true;
                        }
                    }

                    // Has this edge been converted to a physical boundary?
                    const labelList& meFaces = edgeFaces_[meIndex];

                    forAll(meFaces, faceI)
                    {
                        label mfPatch = whichPatch(meFaces[faceI]);

                        if (mfPatch == -1)
                        {
                            continue;
                        }

                        if (getNeighbourProcessor(mfPatch) == -1)
                        {
                            // Store physical patch for conversion
                            newPatch = mfPatch;
                        }
                        else
                        {
                            requireConversion = false;
                        }
                    }

                    // Was an appropriate physical patch found?
                    if (physConvert && !requireConversion)
                    {
                        continue;
                    }

                    if (requireConversion)
                    {
                        edge newEdge = edges_[meIndex];
                        labelList newEdgeFaces = edgeFaces_[meIndex];

                        // Insert the new edge
                        label newEdgeIndex =
                        (
                            insertEdge
                            (
                                newPatch,
                                newEdge,
                                newEdgeFaces
                            )
                        );

                        // Replace faceEdges for all
                        // connected faces.
                        forAll(newEdgeFaces, faceI)
                        {
                            meshOps::replaceLabel
                            (
                                meIndex,
                                newEdgeIndex,
                                faceEdges_[newEdgeFaces[faceI]]
                            );
                        }

                        // Remove the edge
                        removeEdge(meIndex);

                        // Update map
                        map.removeEdge(meIndex);
                        map.addEdge(newEdgeIndex, labelList(1, meIndex));

                        // Replace index and patch
                        meIndex = newEdgeIndex;
                        mePatch = newPatch;

                        // If conversion was to a physical patch,
                        // skip the remaining face mapping steps
                        if (getNeighbourProcessor(newPatch) == -1)
                        {
                            continue;
                        }
                    }

                    // Fetch edges connected to first slave point
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
                                    << " :: " << edges_[meIndex]
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
                        // Check if a coupling existed before
                        if (edgeMap.found(meIndex) && slaveMapPtr)
                        {
                            // Alias for convenience...
                            const changeMap& sMap = *slaveMapPtr;
                            const labelList& reList = sMap.removedEdgeList();

                            // Was the edge removed by the slave?
                            if (findIndex(reList, edgeMap[meIndex]) > -1)
                            {
                                if (debug > 2)
                                {
                                    Pout<< " Found removed edge: "
                                        << edgeMap[meIndex]
                                        << " with meIndex: " << meIndex
                                        << " :: " << edges_[meIndex]
                                        << " on proc: " << procIndices_[pI]
                                        << endl;
                                }

                                // Remove map entry
                                rEdgeMap.erase(edgeMap[meIndex]);
                                edgeMap.erase(meIndex);

                                matchedEdges[edgeI] = true;
                            }
                        }
                    }

                    if ((debug > 4) && !matchedEdges[edgeI])
                    {
                        sMesh.writeVTK
                        (
                            "failedEdgePoints_"
                          + Foam::name(meIndex),
                            labelList(cE), 0, false, true
                        );

                        Pout<< " Failed to match edge: "
                            << meIndex << " :: " << edges_[meIndex]
                            << " using comparison edge: " << cE
                            << " on proc: " << procIndices_[pI]
                            << endl;
                    }
                }
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

    // Set the flag
    topoChangeFlag_ = true;

    // Increment the counter
    status(TOTAL_COLLAPSES)++;

    // Increment the number of modifications
    status(TOTAL_MODIFICATIONS)++;

    // Return a succesful collapse
    map.type() = collapseCase;

    return map;
}


// Remove cell layer above specified patch
const changeMap dynamicTopoFvMesh::removeCellLayer
(
    const label patchID
)
{
    changeMap map;

    labelHashSet edgesToRemove, facesToRemove;
    Map<labelPair> pointsToRemove, edgesToKeep;

    dynamicLabelList patchFaces(patchSizes_[patchID]);
    DynamicList<labelPair> cellsToRemove(patchSizes_[patchID]);
    DynamicList<labelPair> hFacesToRemove(patchSizes_[patchID]);

    for (label faceI = nOldInternalFaces_; faceI < faces_.size(); faceI++)
    {
        label pIndex = whichPatch(faceI);

        if (pIndex != patchID)
        {
            continue;
        }

        // Fetch owner cell
        label cIndex = owner_[faceI];

        // Add face to the list
        patchFaces.append(faceI);

        // Fetch appropriate face / cell
        const face& bFace = faces_[faceI];
        const cell& ownCell = cells_[cIndex];

        // Get the opposing face from the cell
        oppositeFace oFace = ownCell.opposingFace(faceI, faces_);

        if (!oFace.found())
        {
            // Something's wrong here.
            FatalErrorIn
            (
                "const changeMap dynamicTopoFvMesh::removeCellLayer"
                "(const label patchID)"
            )
                << " Face: " << faceI << " :: " << bFace << nl
                << " has no opposing face in cell: "
                << cIndex << " :: " << ownCell << nl
                << abort(FatalError);
        }

        // Fetch cell on the other-side of the opposite face
        label otherCellIndex =
        (
            (owner_[oFace.oppositeIndex()] == cIndex) ?
            neighbour_[oFace.oppositeIndex()] :
            owner_[oFace.oppositeIndex()]
        );

        if (otherCellIndex == -1)
        {
            // Opposite face is on a boundary, and layer
            // removal would be invalid if we continued.
            map.type() = -2;

            return map;
        }

        // Fetch reference to other cell
        const cell& otherCell = cells_[otherCellIndex];

        // Find opposite face on the other cell
        oppositeFace otFace =
        (
            otherCell.opposingFace
            (
                oFace.oppositeIndex(),
                faces_
            )
        );

        if (!otFace.found())
        {
            // Something's wrong here.
            FatalErrorIn
            (
                "const changeMap dynamicTopoFvMesh::removeCellLayer"
                "(const label patchID)"
            )
                << " Face: " << oFace.oppositeIndex()
                << " :: " << oFace << nl
                << " has no opposing face in cell: "
                << otherCellIndex << " :: " << otherCell << nl
                << abort(FatalError);
        }

        // All edges on the boundary face are to be retained
        const labelList& fEdges = faceEdges_[faceI];
        const labelList& ofEdges = faceEdges_[oFace.oppositeIndex()];
        const labelList& otfEdges = faceEdges_[otFace.oppositeIndex()];

        forAll(fEdges, edgeI)
        {
            label eIndex = fEdges[edgeI];

            if (!edgesToKeep.found(eIndex))
            {
                // Find equivalent edge on opposite face
                label oeIndex = -1, oteIndex = -1;
                const edge& bEdge = edges_[eIndex];

                // Build edges for comparison
                label startLoc = bFace.which(bEdge[0]);
                label endLoc = bFace.which(bEdge[1]);

                edge cEdge(oFace[startLoc], oFace[endLoc]);
                edge ctEdge(otFace[startLoc], otFace[endLoc]);

                forAll(ofEdges, edgeJ)
                {
                    const edge& ofEdge = edges_[ofEdges[edgeJ]];

                    if (cEdge == ofEdge)
                    {
                        oeIndex = ofEdges[edgeJ];
                        break;
                    }
                }

                forAll(otfEdges, edgeJ)
                {
                    const edge& otfEdge = edges_[otfEdges[edgeJ]];

                    if (ctEdge == otfEdge)
                    {
                        oteIndex = otfEdges[edgeJ];
                        break;
                    }
                }

                if (oeIndex < 0 || oteIndex < 0)
                {
                    FatalErrorIn
                    (
                        "const changeMap dynamicTopoFvMesh::removeCellLayer"
                        "(const label patchID)"
                    )
                        << " Could not find comparison edge: " << nl
                        << " cEdge: " << cEdge
                        << " oeIndex: " << oeIndex << nl
                        << " ctEdge: " << ctEdge
                        << " oteIndex: " << oteIndex << nl
                        << " for edge: " << bEdge
                        << abort(FatalError);
                }

                // Make entry
                edgesToKeep.insert(eIndex, labelPair(oeIndex, oteIndex));
            }
        }

        // Add information to removal lists
        cellsToRemove.append
        (
            labelPair
            (
                cIndex,
                otherCellIndex
            )
        );

        hFacesToRemove.append
        (
            labelPair
            (
                oFace.oppositeIndex(),
                otFace.oppositeIndex()
            )
        );

        // Mark points for removal
        forAll(oFace, pointI)
        {
            label pIndex = oFace[pointI];

            if (!pointsToRemove.found(pIndex))
            {
                // Make entry
                pointsToRemove.insert
                (
                    pIndex,
                    labelPair(bFace[pointI], otFace[pointI])
                );
            }
        }

        // Loop through cell faces and mark
        // faces / edges for removal
        forAll(ownCell, faceJ)
        {
            label fIndex = ownCell[faceJ];

            if (fIndex == faceI || fIndex == oFace.oppositeIndex())
            {
                continue;
            }

            if (!facesToRemove.found(fIndex))
            {
                facesToRemove.insert(fIndex);
            }

            const labelList& checkEdges = faceEdges_[fIndex];

            forAll(checkEdges, edgeI)
            {
                label eIndex = checkEdges[edgeI];

                if (edgesToKeep.found(eIndex))
                {
                    continue;
                }

                if (!edgesToRemove.found(eIndex))
                {
                    edgesToRemove.insert(eIndex);
                }
            }
        }
    }

    // Correct edgeFaces / faceEdges for retained edges
    forAllConstIter(Map<labelPair>, edgesToKeep, eIter)
    {
        const labelList& rmeFaces = edgeFaces_[eIter().first()];

        forAll(rmeFaces, faceI)
        {
            labelList& fE = faceEdges_[rmeFaces[faceI]];

            bool foundRp = (findIndex(fE, eIter.key()) > -1);
            bool foundRn = (findIndex(fE, eIter().second()) > -1);

            if (foundRp)
            {
                // Size-down edgeFaces for replacement
                meshOps::sizeDownList
                (
                    rmeFaces[faceI],
                    edgeFaces_[eIter.key()]
                );
            }

            if (foundRn)
            {
                // Size-up edgeFaces for replacement
                meshOps::sizeUpList
                (
                    rmeFaces[faceI],
                    edgeFaces_[eIter.key()]
                );

                // Replace edge index
                meshOps::replaceLabel
                (
                    eIter().first(),
                    eIter.key(),
                    fE
                );
            }
        }
    }

    // Remove unwanted faces
    forAllConstIter(labelHashSet, facesToRemove, fIter)
    {
        // Remove the face
        removeFace(fIter.key());

        // Update map
        map.removeFace(fIter.key());
    }

    // Remove unwanted edges
    forAllConstIter(labelHashSet, edgesToRemove, eIter)
    {
        // Remove the edge
        removeEdge(eIter.key());

        // Update map
        map.removeEdge(eIter.key());
    }

    // Remove unwanted points
    forAllConstIter(Map<labelPair>, pointsToRemove, pIter)
    {
        // Update pointEdges information first
        if (is3D())
        {
            const labelList& pEdges = pointEdges_[pIter.key()];

            // Configure edge for comparison
            edge cEdge
            (
                pIter.key(),
                pIter().second()
            );

            label replaceEdge = -1;

            forAll(pEdges, edgeI)
            {
                const edge& checkEdge = edges_[pEdges[edgeI]];

                if (checkEdge == cEdge)
                {
                    replaceEdge = pEdges[edgeI];
                    break;
                }
            }

            if (replaceEdge == -1)
            {
                FatalErrorIn
                (
                    "const changeMap dynamicTopoFvMesh::removeCellLayer"
                    "(const label patchID)"
                )
                    << " Could not find comparison edge: " << nl
                    << " cEdge: " << cEdge
                    << " for point: " << pIter.key() << nl
                    << " pointEdges: " << pEdges
                    << abort(FatalError);
            }

            // Size-up pointEdges
            meshOps::sizeUpList
            (
                replaceEdge,
                pointEdges_[pIter().first()]
            );
        }

        // Remove the point
        removePoint(pIter.key());

        // Update map
        map.removePoint(pIter.key());
    }

    // Track modified faces for mapping
    labelHashSet modifiedFaces;

    // Remove all cells
    forAll(cellsToRemove, indexI)
    {
        // Replace face label on the other cell
        meshOps::replaceLabel
        (
            hFacesToRemove[indexI].first(),
            patchFaces[indexI],
            cells_[cellsToRemove[indexI].second()]
        );

        // Set owner information
        owner_[patchFaces[indexI]] = cellsToRemove[indexI].second();

        // Replace points on faces / edges
        const cell& otherCell = cells_[cellsToRemove[indexI].second()];

        forAll(otherCell, faceI)
        {
            face& faceToCheck = faces_[otherCell[faceI]];

            bool modified = false;

            forAll(faceToCheck, pointI)
            {
                if (pointsToRemove.found(faceToCheck[pointI]))
                {
                    faceToCheck[pointI] =
                    (
                        pointsToRemove[faceToCheck[pointI]].first()
                    );
                }

                modified = true;
            }

            // Add face to list
            if (modified && !modifiedFaces.found(otherCell[faceI]))
            {
                modifiedFaces.insert(otherCell[faceI]);
            }

            const labelList& fEdges = faceEdges_[otherCell[faceI]];

            forAll(fEdges, edgeI)
            {
                edge& edgeToCheck = edges_[fEdges[edgeI]];

                if (pointsToRemove.found(edgeToCheck[0]))
                {
                    edgeToCheck[0] =
                    (
                        pointsToRemove[edgeToCheck[0]].first()
                    );
                }

                if (pointsToRemove.found(edgeToCheck[1]))
                {
                    edgeToCheck[1] =
                    (
                        pointsToRemove[edgeToCheck[1]].first()
                    );
                }
            }
        }

        // Remove horizontal interior face
        removeFace(hFacesToRemove[indexI].first());

        // Update map
        map.removeFace(hFacesToRemove[indexI].first());

        // Remove the cell
        removeCell(cellsToRemove[indexI].first());

        // Update map
        map.removeCell(cellsToRemove[indexI].first());
    }

    // Now that all old / new cells possess correct connectivity,
    // compute mapping information.
    forAll(cellsToRemove, indexI)
    {
        // Set mapping for the modified cell
        setCellMapping
        (
            cellsToRemove[indexI].second(),
            labelList(cellsToRemove[indexI])
        );
    }

    // Set face mapping information for modified faces
    forAllConstIter(labelHashSet, modifiedFaces, fIter)
    {
        // Decide between default / weighted mapping
        // based on boundary information
        label fPatch = whichPatch(fIter.key());

        if (fPatch == -1)
        {
            // Default mapping for interior faces
            setFaceMapping(fIter.key());
        }
        else
        {
            // Map boundary face from itself
            setFaceMapping(fIter.key(), labelList(1, fIter.key()));
        }
    }

    // Set the flag
    topoChangeFlag_ = true;

    // Specify that the operation was successful
    map.type() = 1;

    // Return the changeMap
    return map;
}


// Merge a set of boundary faces into internal
const changeMap dynamicTopoFvMesh::mergeBoundaryFaces
(
    const labelList& mergeFaces
)
{
    changeMap map;

    if (debug > 2)
    {
        Pout << "Merging faces: " << mergeFaces << endl;
    }

    List<labelPair> mergePairs(0);

    forAll(mergeFaces, faceI)
    {
        const face& fI = faces_[mergeFaces[faceI]];

        for (label faceJ = faceI + 1; faceJ < mergeFaces.size(); faceJ++)
        {
            const face& fJ = faces_[mergeFaces[faceJ]];

            bool merge = false;

            if (fI.size() == fJ.size() && fI.size() == 3)
            {
                if
                (
                    triFace::compare
                    (
                        triFace(fI[0], fI[1], fI[2]),
                        triFace(fJ[0], fJ[1], fJ[2])
                    )
                )
                {
                    merge = true;
                }
            }
            else
            if (face::compare(fI, fJ))
            {
                merge = true;
            }

            if (merge)
            {
                meshOps::sizeUpList
                (
                    labelPair(mergeFaces[faceI], mergeFaces[faceJ]),
                    mergePairs
                );

                break;
            }
        }
    }

    if (debug > 2)
    {
        Pout<< " mergePairs: " << mergePairs << endl;
    }

    dynamicLabelList checkEdges(10);

    forAll(mergePairs, pairI)
    {
        label firstFace = mergePairs[pairI].first();
        label secondFace = mergePairs[pairI].second();

        // Obtain owners for both faces, and compare their labels
        label newOwner = -1, newNeighbour = -1;
        label removedFace = -1, retainedFace = -1;

        if (owner_[firstFace] < owner_[secondFace])
        {
            // Retain the first face
            removedFace = secondFace;
            retainedFace = firstFace;

            newOwner = owner_[firstFace];
            newNeighbour = owner_[secondFace];
        }
        else
        {
            // Retain the second face
            removedFace = firstFace;
            retainedFace = secondFace;

            newOwner = owner_[secondFace];
            newNeighbour = owner_[firstFace];
        }

        // Insert a new interior face
        label newFaceIndex =
        (
            insertFace
            (
                -1,
                face(faces_[retainedFace]),
                newOwner,
                newNeighbour,
                labelList(faceEdges_[retainedFace])
            )
        );

        // Update map
        map.addFace(newFaceIndex);

        // Replace cell with the new face label
        meshOps::replaceLabel
        (
            removedFace,
            newFaceIndex,
            cells_[owner_[removedFace]]
        );

        meshOps::replaceLabel
        (
            retainedFace,
            newFaceIndex,
            cells_[owner_[retainedFace]]
        );

        const labelList& fEdges = faceEdges_[newFaceIndex];
        const labelList& rfEdges = faceEdges_[removedFace];

        // Check for common edges on the removed face
        forAll(rfEdges, edgeI)
        {
            label reIndex = rfEdges[edgeI];

            if (findIndex(fEdges, reIndex) == -1)
            {
                // Find the equivalent edge
                const edge& rEdge = edges_[reIndex];
                const labelList& reFaces = edgeFaces_[reIndex];

                label keIndex = -1;

                forAll(fEdges, edgeJ)
                {
                    if (edges_[fEdges[edgeJ]] == rEdge)
                    {
                        keIndex = fEdges[edgeJ];
                        break;
                    }
                }

                if (keIndex == -1)
                {
                    Pout<< " Could not find edge for: "
                        << reIndex << " :: " << rEdge
                        << abort(FatalError);
                }

                // Add edgeFaces from this edge
                forAll(reFaces, faceI)
                {
                    if (reFaces[faceI] == removedFace)
                    {
                        continue;
                    }

                    meshOps::sizeUpList
                    (
                        reFaces[faceI],
                        edgeFaces_[keIndex]
                    );

                    meshOps::replaceLabel
                    (
                        reIndex,
                        keIndex,
                        faceEdges_[reFaces[faceI]]
                    );
                }

                // Remove the old edge
                removeEdge(reIndex);

                // Update map
                map.removeEdge(reIndex);
            }
        }

        // Replace edgeFaces with the new face label
        forAll(fEdges, edgeI)
        {
            label eIndex = fEdges[edgeI];

            if (findIndex(edgeFaces_[eIndex], removedFace) > -1)
            {
                meshOps::sizeDownList
                (
                    removedFace,
                    edgeFaces_[eIndex]
                );
            }

            if (findIndex(edgeFaces_[eIndex], retainedFace) > -1)
            {
                meshOps::sizeDownList
                (
                    retainedFace,
                    edgeFaces_[eIndex]
                );
            }

            // Size-up list with the new face index
            meshOps::sizeUpList
            (
                newFaceIndex,
                edgeFaces_[eIndex]
            );

            if (findIndex(checkEdges, eIndex) == -1)
            {
                checkEdges.append(eIndex);
            }
        }

        // Remove the boundary faces
        removeFace(removedFace);
        removeFace(retainedFace);

        // Update map
        map.removeFace(removedFace);
        map.removeFace(retainedFace);
    }

    forAll(checkEdges, edgeI)
    {
        bool allInterior = true;
        label eIndex = checkEdges[edgeI];

        const labelList& eFaces = edgeFaces_[eIndex];

        forAll(eFaces, faceI)
        {
            if (whichPatch(eFaces[faceI]) > -1)
            {
                allInterior = false;
                break;
            }
        }

        if (allInterior)
        {
            // This edge needs to be converted to an interior one
            label newEdgeIndex =
            (
                insertEdge
                (
                    -1,
                    edge(edges_[eIndex]),
                    labelList(edgeFaces_[eIndex])
                )
            );

            // Update map
            map.addEdge(newEdgeIndex, labelList(1, eIndex));

            // Update faceEdges information for all connected faces
            const labelList& neFaces = edgeFaces_[newEdgeIndex];

            forAll(neFaces, faceI)
            {
                meshOps::replaceLabel
                (
                    eIndex,
                    newEdgeIndex,
                    faceEdges_[neFaces[faceI]]
                );
            }

            // Remove the old boundary edge
            removeEdge(eIndex);

            // Update map
            map.removeEdge(eIndex);

            // Replace index
            checkEdges[edgeI] = newEdgeIndex;
        }
    }

    if (debug > 2)
    {
        Pout<< "Merge complete." << nl << endl;
    }

    // Return a succesful merge
    map.type() = 1;

    return map;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
