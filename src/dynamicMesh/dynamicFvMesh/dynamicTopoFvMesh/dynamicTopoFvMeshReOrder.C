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

#include "objectMap.H"
#include "dynamicTopoFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Reorder points after a topology change
void dynamicTopoFvMesh::reOrderPoints
(
    pointField& points,
    pointField& preMotionPoints,
    labelListList& pointZoneMap,
    bool threaded
)
{
    // *** Point renumbering *** //
    // If points were deleted during topology change, the numerical order ceases
    // to be continuous. Loop through all points and renumber sequentially.

    if (debug)
    {
        if (threaded)
        {
            Info << "Thread: " << self() << ": ";
        }

        Info << "ReOrdering points..." << flush;
    }

    // Allocate for the mapping information
    pointMap_.setSize(nPoints_, -1);

    label pointInOrder = 0;

    addedPointRenumbering_.clear();

    for (label pointI = 0; pointI < nOldPoints_; pointI++)
    {
        // Check if this is a deleted point
        if (reversePointMap_[pointI] == -1)
        {
            continue;
        }

        // Update the point info
        points[pointInOrder] = points_[pointI];
        preMotionPoints[pointInOrder] = oldPoints_[pointI];

        // Update maps
        pointMap_[pointInOrder]  = pointI;
        reversePointMap_[pointI] = pointInOrder;

        // Update the counter
        pointInOrder++;
    }

    for (label pointI = nOldPoints_; pointI < points_.size(); pointI++)
    {
        // Was this point removed after addition?
        if (deletedPoints_.found(pointI))
        {
            continue;
        }

        // Update the point info
        points[pointInOrder] = points_[pointI];
        preMotionPoints[pointInOrder] = oldPoints_[pointI];

        // Put inserted points in a seperate hashSet
        addedPointRenumbering_.insert(pointI, pointInOrder);

        // Update the counter
        pointInOrder++;
    }

    // Now that we're done preparing the point maps, unlock the point mutex
    if (threaded)
    {
        entityMutex_[0].unlock();
    }

    // Final check to ensure everything went okay
    if (pointInOrder != nPoints_)
    {
        FatalErrorIn("dynamicTopoFvMesh::reOrderPoints()") << nl
            << " Algorithm did not visit every point in the mesh."
            << " Something's messed up." << nl
            << abort(FatalError);
    }

    // Wipe out pointsFromPoints, since this is not
    // really used in this context for mapping.
    pointsFromPoints_.clear();

    // Renumber all maps.
    forAll(pointsFromPoints_, indexI)
    {
        objectMap& thisMap = pointsFromPoints_[indexI];

        if (thisMap.index() < nOldPoints_)
        {
            thisMap.index() = reversePointMap_[thisMap.index()];
        }
        else
        {
            thisMap.index() = addedPointRenumbering_[thisMap.index()];
        }
    }

    // Prepare the pointZoneMap
    pointZoneMesh& pointZones = polyMesh::pointZones();

    labelListList newPointZoneAddr(pointZones.size());
    labelList nPointsInZone(pointZones.size(), 0);

    // Prepare zone maps.
    forAll(pointZones, pzI)
    {
        // Get the list of old points
        const labelList& oldAddr = pointZones[pzI];

        label& curNPoints = nPointsInZone[pzI];

        // First count the actual number of points in each zone.
        forAll(oldAddr, pointI)
        {
            // Was this zone point deleted? Don't count it.
            if (reversePointMap_[oldAddr[pointI]] != -1)
            {
                curNPoints++;
            }
        }

        // Check for added points as well
        forAllIter(Map<label>, addedPointZones_, pIter)
        {
            if (pIter() == pzI)
            {
                curNPoints++;
            }
        }

        label pIndex = 0;
        labelList& newAddr = newPointZoneAddr[pzI];

        // Set the sizes first
        newAddr.setSize(curNPoints);
        pointZoneMap[pzI].setSize(curNPoints, -1);

        // Add existing zone points which have been renumbered.
        forAll(oldAddr, pointI)
        {
            if (reversePointMap_[oldAddr[pointI]] != -1)
            {
                pointZoneMap[pzI][pIndex] = oldAddr[pointI];
                newAddr[pIndex] = reversePointMap_[oldAddr[pointI]];
                pIndex++;
            }
        }

        // Next, add the newly added zone points.
        forAllIter(Map<label>, addedPointZones_, pIter)
        {
            if (pIter() == pzI)
            {
                newAddr[pIndex++] = addedPointRenumbering_[pIter.key()];
            }
        }

        // Finally, assign addressing to this zone.
        pointZones[pzI] = newPointZoneAddr[pzI];
    }

    // Reset all zones
    pointZones.updateMesh();

    // Clear local point copies
    points_.clear();
    oldPoints_.clear();

    // Set values
    points_ = points;
    oldPoints_ = preMotionPoints;

    if (debug)
    {
        Info << "Done." << endl;
    }
}


// Static equivalent for multi-threading
void dynamicTopoFvMesh::reOrderPointsThread
(
    void *argument
)
{
    // Recast the argument
    meshHandler *thread = static_cast<meshHandler*>(argument);

    dynamicTopoFvMesh& mesh = thread->reference();

    // Lock the point mutex first
    mesh.entityMutex(0).lock();

    // Signal the calling thread
    thread->sendSignal(meshHandler::START);

    // Recast the pointers for the reOrderPoints argument
    pointField& points =
    (
        *(static_cast<pointField*>(thread->operator()(0)))
    );

    pointField& preMotionPoints =
    (
        *(static_cast<pointField*>(thread->operator()(1)))
    );

    labelListList& pointZoneMap =
    (
        *(static_cast<labelListList*>(thread->operator()(2)))
    );

    // Reorder the points
    mesh.reOrderPoints(points, preMotionPoints, pointZoneMap, true);
}


// Reorder edges after a topology change
void dynamicTopoFvMesh::reOrderEdges
(
    edgeList& edges,
    labelListList& edgeFaces,
    labelListList& faceEdges,
    bool threaded
)
{
    // *** Edge renumbering *** //
    // If edges were deleted during topology change, the numerical order ceases
    // to be continuous. Edges are added to respective internal/boundary patches

    if (debug)
    {
        if (threaded)
        {
            Info << "Thread: " << self() << ": ";
        }

        Info << "ReOrdering edges..." << flush;
    }

    // Allocate for mapping information
    edgeMap_.setSize(nEdges_, -1);

    label edgeInOrder = 0, allEdges = edges_.size();
    edgeList oldEdges(allEdges);
    labelListList oldEdgeFaces(allEdges);
    labelListList oldEdgePoints(allEdges);

    addedEdgeRenumbering_.clear();
    Map<label> addedEdgeReverseRenumbering;

    // Transfer old edge-based lists, and clear them
    forAll(edges_, edgeI)
    {
        oldEdges[edgeI] = edges_[edgeI];
        oldEdgeFaces[edgeI].transfer(edgeFaces_[edgeI]);
    }

    edges_.setSize(nEdges_); edgeFaces_.setSize(nEdges_);

    if (!twoDMesh_)
    {
        forAll(edgePoints_, edgeI)
        {
            oldEdgePoints[edgeI].transfer(edgePoints_[edgeI]);
        }

        edgePoints_.setSize(nEdges_);
    }

    // Keep track of inserted boundary edge indices
    labelList boundaryPatchIndices(edgePatchStarts_);

    // Loop through all edges and add internal ones first
    forAll(oldEdges, edgeI)
    {
        // Ensure that we're adding valid edges
        if (oldEdgeFaces[edgeI].empty())
        {
            continue;
        }

        // Determine which patch this edge belongs to
        label patch = whichEdgePatch(edgeI);

        // Update maps for boundary edges. Edge insertion for
        // boundaries will be done after internal edges.
        if (patch >= 0)
        {
            label bEdgeIndex = boundaryPatchIndices[patch]++;

            // Update the maps
            if (edgeI < nOldEdges_)
            {
                edgeMap_[bEdgeIndex] = edgeI;
                reverseEdgeMap_[edgeI] = bEdgeIndex;
            }
            else
            {
                addedEdgeRenumbering_.insert(edgeI, bEdgeIndex);
                addedEdgeReverseRenumbering.insert(bEdgeIndex, edgeI);

                edgeMap_[bEdgeIndex] = -1;
            }
        }
        else
        {
            // Obtain references
            edge& thisEdge = oldEdges[edgeI];
            labelList& thisEF = oldEdgeFaces[edgeI];

            // Renumber internal edges and add normally.
            if (edgeI < nOldEdges_)
            {
                edgeMap_[edgeInOrder] = edgeI;
                reverseEdgeMap_[edgeI] = edgeInOrder;
            }
            else
            {
                addedEdgeRenumbering_.insert(edgeI, edgeInOrder);
            }

            // Insert entities into local lists...
            edges_[edgeInOrder] = thisEdge;
            edgeFaces_[edgeInOrder] = thisEF;

            // Insert entities into mesh-reset lists
            edges[edgeInOrder] = thisEdge;
            edgeFaces[edgeInOrder].transfer(thisEF);

            if (!twoDMesh_)
            {
                edgePoints_[edgeInOrder].transfer(oldEdgePoints[edgeI]);
            }

            edgeInOrder++;
        }
    }

    // All internal edges have been inserted. Now insert boundary edges.
    label oldIndex;

    for(label i = nInternalEdges_; i < nEdges_; i++)
    {
        if (edgeMap_[i] == -1)
        {
            // This boundary edge was added during the topology change
            oldIndex = addedEdgeReverseRenumbering[i];
        }
        else
        {
            oldIndex = edgeMap_[i];
        }

        // Insert entities into local Lists...
        edges_[edgeInOrder] = oldEdges[oldIndex];
        edgeFaces_[edgeInOrder] = oldEdgeFaces[oldIndex];

        // Insert entities into mesh-reset lists
        edges[edgeInOrder] = oldEdges[oldIndex];
        edgeFaces[edgeInOrder].transfer(oldEdgeFaces[oldIndex]);

        if (!twoDMesh_)
        {
            edgePoints_[edgeInOrder].transfer(oldEdgePoints[oldIndex]);
        }

        edgeInOrder++;
    }

    // Now that we're done with edges, unlock it
    if (threaded)
    {
        entityMutex_[1].unlock();
    }

    // Final check to ensure everything went okay
    if (edgeInOrder != nEdges_)
    {
        FatalErrorIn("dynamicTopoFvMesh::reOrderEdges()") << nl
                << " Algorithm did not visit every edge in the mesh."
                << " Something's messed up." << nl
                << abort(FatalError);
    }

    // Renumber all edges / edgePoints with updated point information
    label pIndex = -1;

    if (threaded)
    {
        entityMutex_[0].lock();
    }

    forAll(edges_, edgeI)
    {
        // Obtain references
        edge& thisEdge = edges_[edgeI];
        edge& thisREdge = edges[edgeI];

        // Renumber edges
        if (thisEdge[0] < nOldPoints_)
        {
            pIndex = reversePointMap_[thisEdge[0]];
        }
        else
        {
            pIndex = addedPointRenumbering_[thisEdge[0]];
        }

        thisEdge[0] = pIndex;
        thisREdge[0] = pIndex;

        if (thisEdge[1] < nOldPoints_)
        {
            pIndex = reversePointMap_[thisEdge[1]];
        }
        else
        {
            pIndex = addedPointRenumbering_[thisEdge[1]];
        }

        thisEdge[1] = pIndex;
        thisREdge[1] = pIndex;

        // Renumber edgePoints
        if (!twoDMesh_)
        {
            labelList& ePoints = edgePoints_[edgeI];

            forAll(ePoints, pointI)
            {
                if (ePoints[pointI] < nOldPoints_)
                {
                    ePoints[pointI] = reversePointMap_[ePoints[pointI]];
                }
                else
                {
                    ePoints[pointI] = addedPointRenumbering_[ePoints[pointI]];
                }
            }
        }
    }

    if (threaded)
    {
        entityMutex_[0].unlock();
    }

    // Renumber all faceEdges
    label eIndex = -1;

    if (threaded)
    {
        entityMutex_[2].lock();
    }

    forAll(faceEdges_, faceI)
    {
        // Obtain references
        labelList& fEdges = faceEdges_[faceI];
        labelList& rfEdges = faceEdges[faceI];

        forAll(fEdges, edgeI)
        {
            if (fEdges[edgeI] < nOldEdges_)
            {
                eIndex = reverseEdgeMap_[fEdges[edgeI]];
            }
            else
            {
                eIndex = addedEdgeRenumbering_[fEdges[edgeI]];
            }

            fEdges[edgeI] = eIndex;
            rfEdges[edgeI] = eIndex;
        }
    }

    // Renumber all edgeFaces
    label fIndex = -1;

    forAll(edgeFaces_, edgeI)
    {
        // Obtain references
        labelList& eFaces = edgeFaces_[edgeI];
        labelList& reFaces = edgeFaces[edgeI];

        // Renumber edgeFaces
        forAll(eFaces, faceI)
        {
            if (eFaces[faceI] < nOldFaces_)
            {
                fIndex = reverseFaceMap_[eFaces[faceI]];
            }
            else
            {
                fIndex = addedFaceRenumbering_[eFaces[faceI]];
            }

            eFaces[faceI] = fIndex;
            reFaces[faceI] = fIndex;
        }
    }

    if (threaded)
    {
        entityMutex_[2].unlock();
    }

    if (!twoDMesh_)
    {
        // Invert edges to obtain pointEdges
        pointEdges_ = invertManyToMany<edge, labelList>(nPoints_, edges_);
    }

    if (debug)
    {
        Info << "Done." << endl;
    }
}


// Static equivalent for multi-threading
void dynamicTopoFvMesh::reOrderEdgesThread
(
    void *argument
)
{
    // Recast the argument
    meshHandler *thread = static_cast<meshHandler*>(argument);

    dynamicTopoFvMesh& mesh = thread->reference();

    // Lock the edge mutex first
    mesh.entityMutex(1).lock();

    // Signal the calling thread
    thread->sendSignal(meshHandler::START);

    // Recast the pointers for the argument
    edgeList& edges =
    (
        *(static_cast<edgeList*>(thread->operator()(0)))
    );

    labelListList& edgeFaces =
    (
        *(static_cast<labelListList*>(thread->operator()(1)))
    );

    labelListList& faceEdges =
    (
        *(static_cast<labelListList*>(thread->operator()(2)))
    );

    // Reorder the edges
    mesh.reOrderEdges(edges, edgeFaces, faceEdges, true);

    // Signal the calling thread
    thread->sendSignal(meshHandler::STOP);
}


// Reorder faces in upper-triangular order after a topology change
void dynamicTopoFvMesh::reOrderFaces
(
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    labelListList& faceEdges,
    labelListList& faceZoneFaceMap,
    bool threaded
)
{
    // *** Face renumbering *** //
    // Faces have to be renumbered if any were added/deleted/modified
    // Boundary faces are added to respective patches.
    // Internal faces, however, have to be added in upper-triangular ordering;
    // i.e., in the increasing order of neighbours

    if (debug)
    {
        if (threaded)
        {
            Info << "Thread: " << self() << ": ";
        }

        Info << "ReOrdering faces..." << flush;
    }

    // Allocate for mapping information
    faceMap_.setSize(nFaces_, -1);

    label faceInOrder = 0, allFaces = faces_.size();
    faceList oldFaces(allFaces);
    labelList oldOwner(allFaces), oldNeighbour(allFaces), visited(allFaces,0);
    labelListList oldFaceEdges(allFaces);

    addedFaceRenumbering_.clear();
    Map<label> addedFaceReverseRenumbering;

    // Make a copy of the old face-based lists, and clear them
    forAll(faces_, faceI)
    {
        oldFaces[faceI].transfer(faces_[faceI]);
        oldOwner[faceI] = owner_[faceI];
        oldNeighbour[faceI] = neighbour_[faceI];
        oldFaceEdges[faceI].transfer(faceEdges_[faceI]);
    }

    // Renumber all faces with updated point information
    label pIndex = -1;

    if (threaded)
    {
        entityMutex_[0].lock();
    }

    forAll(oldFaces, faceI)
    {
        face& thisFace = oldFaces[faceI];

        forAll(thisFace, pointI)
        {
            if (thisFace[pointI] < nOldPoints_)
            {
                pIndex = reversePointMap_[thisFace[pointI]];
            }
            else
            {
                pIndex = addedPointRenumbering_[thisFace[pointI]];
            }

            thisFace[pointI] = pIndex;
        }
    }

    if (threaded)
    {
        entityMutex_[0].unlock();
    }

    // Wait for the cell mutex to become available
    if (threaded)
    {
        entityMutex_[3].lock();
    }

    faces_.setSize(nFaces_);
    owner_.setSize(nFaces_);
    neighbour_.setSize(nFaces_);
    faceEdges_.setSize(nFaces_);

    // Mark the internal faces with -2 so that they are inserted first
    forAll(cells_, cellI)
    {
        const cell& curFaces = cells_[cellI];

        forAll(curFaces, faceI)
        {
            visited[curFaces[faceI]]--;
        }
    }

    // Keep track of inserted boundary face indices
    labelList boundaryPatchIndices(patchStarts_);

    // Handle boundaries first. If any coupled interfaces need to be
    // updated, they can be reshuffled after interior faces are done.
    // Update maps for boundaries now.
    for (label faceI = nOldInternalFaces_; faceI < allFaces; faceI++)
    {
        if (visited[faceI] == -1)
        {
            label patchID = whichPatch(faceI);
            label bFaceIndex = boundaryPatchIndices[patchID]++;

            // Update the maps
            if (faceI < nOldFaces_)
            {
                faceMap_[bFaceIndex] = faceI;
                reverseFaceMap_[faceI] = bFaceIndex;
            }
            else
            {
                addedFaceRenumbering_.insert(faceI, bFaceIndex);

                addedFaceReverseRenumbering.insert(bFaceIndex, faceI);

                faceMap_[bFaceIndex] = -1;
            }

            // Mark this face as visited
            visited[faceI] = 0;
        }
    }

    // Upper-triangular ordering of internal faces:

    // Insertion cannot be done in one go as the faces need to be
    // added into the list in the increasing order of neighbour
    // cells.  Therefore, all neighbours will be detected first
    // and then added in the correct order.
    forAll(cells_, cellI)
    {
        // Record the neighbour cell
        const cell& curFaces = cells_[cellI];

        labelList neiCells(curFaces.size(), -1);

        label nNeighbours = 0;

        forAll(curFaces, faceI)
        {
            if (visited[curFaces[faceI]] == -2)
            {
                // Face is internal and gets reordered
                label own =
                (
                    oldOwner[curFaces[faceI]] < nOldCells_
                  ? reverseCellMap_[oldOwner[curFaces[faceI]]]
                  : addedCellRenumbering_[oldOwner[curFaces[faceI]]]
                );

                label nei =
                (
                    oldNeighbour[curFaces[faceI]] < nOldCells_
                  ? reverseCellMap_[oldNeighbour[curFaces[faceI]]]
                  : addedCellRenumbering_[oldNeighbour[curFaces[faceI]]]
                );

                label smallerIndex = own < nei ? own : nei;
                label largerIndex  = own > nei ? own : nei;

                if (cellI == smallerIndex)
                {
                    neiCells[faceI] = largerIndex;
                    nNeighbours++;
                }
            }
        }

        // Add internal faces in the increasing order of neighbours
        for (label neiSearch = 0; neiSearch < nNeighbours; neiSearch++)
        {
            // Find the lowest neighbour which is still valid
            label nextNei = -1;
            label minNei = nCells_;

            forAll(neiCells, ncI)
            {
                if (neiCells[ncI] > -1 && neiCells[ncI] < minNei)
                {
                    nextNei = ncI;
                    minNei = neiCells[ncI];
                }
            }

            if (nextNei > -1)
            {
                // Face is internal and gets reordered
                if (curFaces[nextNei] < nOldFaces_)
                {
                    faceMap_[faceInOrder] = curFaces[nextNei];
                    reverseFaceMap_[curFaces[nextNei]] = faceInOrder;
                }
                else
                {
                    addedFaceRenumbering_.insert
                    (
                        curFaces[nextNei],
                        faceInOrder
                    );
                }

                // Renumber owner/neighbour
                label oldOwn = oldOwner[curFaces[nextNei]];
                label oldNei = oldNeighbour[curFaces[nextNei]];

                label ownerRenumber =
                (
                    oldOwn < nOldCells_
                  ? reverseCellMap_[oldOwn] : addedCellRenumbering_[oldOwn]
                );

                label neighbourRenumber =
                (
                    oldNei < nOldCells_
                  ? reverseCellMap_[oldNei] : addedCellRenumbering_[oldNei]
                );

                // Cell-reordering may cause flipped faces.. Correct them.
                face& faceRenumber = oldFaces[curFaces[nextNei]];

                if (neighbourRenumber < ownerRenumber)
                {
                    faceRenumber = faceRenumber.reverseFace();

                    if (!bandWidthReduction_)
                    {
                        FatalErrorIn("dynamicTopoFvMesh::reOrderFaces()")
                            << nl
                            << " Found an improperly ordered face."
                            << " Something's messed up." << nl
                            << abort(FatalError);
                    }

                    setFlip(curFaces[nextNei]);
                }

                // Insert entities into local lists...
                faces_[faceInOrder] = faceRenumber;
                owner_[faceInOrder] = cellI;
                neighbour_[faceInOrder] = minNei;
                faceEdges_[faceInOrder] =
                (
                    oldFaceEdges[curFaces[nextNei]]
                );

                // Insert entities into mesh-reset lists
                faces[faceInOrder].transfer(faceRenumber);
                owner[faceInOrder] = cellI;
                neighbour[faceInOrder] = minNei;
                faceEdges[faceInOrder].transfer
                (
                    oldFaceEdges[curFaces[nextNei]]
                );

                // Stop the neighbour from being used again
                neiCells[nextNei] = -1;

                // Mark this face as visited
                visited[curFaces[nextNei]] = 0;

                faceInOrder++;
            }
            else
            {
                FatalErrorIn("dynamicTopoFvMesh::reOrderFaces()") << nl
                    << "Error in internal face insertion" << nl
                    << abort(FatalError);
            }
        }
    }

    // All internal faces have been inserted. Now insert boundary faces.
    label oldIndex;

    for (label i = nInternalFaces_; i < nFaces_; i++)
    {
        if (faceMap_[i] == -1)
        {
            // This boundary face was added during the topology change
            oldIndex = addedFaceReverseRenumbering[i];
        }
        else
        {
            oldIndex = faceMap_[i];
        }

        // Renumber owner/neighbour
        label ownerRenumber =
        (
            oldOwner[oldIndex] < nOldCells_
          ? reverseCellMap_[oldOwner[oldIndex]]
          : addedCellRenumbering_[oldOwner[oldIndex]]
        );

        // Insert entities into local listsLists...
        faces_[faceInOrder] = oldFaces[oldIndex];
        owner_[faceInOrder] = ownerRenumber;
        neighbour_[faceInOrder] = -1;
        faceEdges_[faceInOrder] = oldFaceEdges[oldIndex];

        // Insert entities into mesh-reset lists
        // NOTE: From OF-1.5 onwards, neighbour array
        //       does not store -1 for boundary faces
        faces[faceInOrder].transfer(oldFaces[oldIndex]);
        owner[faceInOrder] = ownerRenumber;
        faceEdges[faceInOrder].transfer(oldFaceEdges[oldIndex]);

        faceInOrder++;
    }

    // Now that we're done with faces, unlock it
    if (threaded)
    {
        entityMutex_[2].unlock();
    }

    // Renumber all maps.
    forAll(facesFromPoints_, indexI)
    {
        objectMap& thisMap = facesFromPoints_[indexI];

        if (thisMap.index() < nOldFaces_)
        {
            thisMap.index() = reverseFaceMap_[thisMap.index()];
        }
        else
        {
            thisMap.index() = addedFaceRenumbering_[thisMap.index()];
        }
    }

    forAll(facesFromEdges_, indexI)
    {
        objectMap& thisMap = facesFromEdges_[indexI];

        if (thisMap.index() < nOldFaces_)
        {
            thisMap.index() = reverseFaceMap_[thisMap.index()];
        }
        else
        {
            thisMap.index() = addedFaceRenumbering_[thisMap.index()];
        }
    }

    forAll(facesFromFaces_, indexI)
    {
        objectMap& thisMap = facesFromFaces_[indexI];

        if (thisMap.index() < nOldFaces_)
        {
            thisMap.index() = reverseFaceMap_[thisMap.index()];
        }
        else
        {
            thisMap.index() = addedFaceRenumbering_[thisMap.index()];
        }
    }

    // Renumber all flipFaces
    labelHashSet flipFaces;

    forAllIter(labelHashSet, flipFaces_, fIter)
    {
        if (fIter.key() < nOldFaces_)
        {
            flipFaces.insert(reverseFaceMap_[fIter.key()]);
        }
        else
        {
            // Added faces cannot be flipped.
            FatalErrorIn("dynamicTopoFvMesh::reOrderFaces()") << nl
                << " Face: " << fIter.key()
                << " is new, and shouldn't be flipped." << nl
                << " nOldFaces: " << nOldFaces_
                << abort(FatalError);
        }
    }

    flipFaces_.transfer(flipFaces);

    // Renumber all cells with updated face information
    forAll(cells_, cellI)
    {
        cell& cellFaces = cells_[cellI];

        forAll(cellFaces, faceI)
        {
            if (cellFaces[faceI] < nOldFaces_)
            {
                cellFaces[faceI] = reverseFaceMap_[cellFaces[faceI]];
            }
            else
            {
                cellFaces[faceI] = addedFaceRenumbering_[cellFaces[faceI]];
            }
        }
    }

    // Now that we're done with cells, unlock it
    if (threaded)
    {
        entityMutex_[3].unlock();
    }

    // Final check to ensure everything went okay
    if (debug > 1)
    {
        if (sum(visited) != 0)
        {
            FatalErrorIn("dynamicTopoFvMesh::reOrderFaces()") << nl
                << " Algorithm did not visit every face in the mesh."
                << " Something's messed up." << nl
                << abort(FatalError);
        }
    }

    // Prepare the faceZoneMap
    faceZoneMesh& faceZones = polyMesh::faceZones();

    labelListList newFaceZoneAddr(faceZones.size());
    boolListList faceZoneFlipMap(faceZones.size());
    labelList nFacesInZone(faceZones.size(), 0);

    // Prepare zone maps.
    forAll(faceZones, fzI)
    {
        // Get the list of old faces
        const labelList& oldAddr = faceZones[fzI];

        label& curNFaces = nFacesInZone[fzI];

        // First count the actual number of faces in each zone.
        forAll(oldAddr, faceI)
        {
            // Was this zone face deleted? Don't count it.
            if (reverseFaceMap_[oldAddr[faceI]] != -1)
            {
                curNFaces++;
            }
        }

        // Check for added faces as well
        forAllIter(Map<label>, addedFaceZones_, fIter)
        {
            if (fIter() == fzI)
            {
                curNFaces++;
            }
        }

        label fIndex = 0;
        labelList& newAddr = newFaceZoneAddr[fzI];

        // Set the sizes first
        newAddr.setSize(curNFaces);
        faceZoneFaceMap[fzI].setSize(curNFaces, -1);
        faceZoneFlipMap[fzI].setSize(curNFaces, false);

        // Add existing zone faces which have been renumbered.
        forAll(oldAddr, faceI)
        {
            if (reverseFaceMap_[oldAddr[faceI]] != -1)
            {
                faceZoneFaceMap[fzI][fIndex] = oldAddr[faceI];
                newAddr[fIndex] = reverseFaceMap_[oldAddr[faceI]];
                fIndex++;
            }
        }

        // Next, add the newly added zone faces.
        forAllIter(Map<label>, addedFaceZones_, fIter)
        {
            if (fIter() == fzI)
            {
                newAddr[fIndex++] = addedFaceRenumbering_[fIter.key()];
            }
        }

        // Finally, reset addressing for this zone.
        faceZones[fzI].resetAddressing
        (
            newFaceZoneAddr[fzI],
            faceZoneFlipMap[fzI]
        );
    }

    // Reset all zones
    faceZones.updateMesh();

    if (debug)
    {
        Info << "Done." << endl;
    }
}


// Static equivalent for multi-threading
void dynamicTopoFvMesh::reOrderFacesThread
(
    void *argument
)
{
    // Recast the argument
    meshHandler *thread = static_cast<meshHandler*>(argument);

    dynamicTopoFvMesh& mesh = thread->reference();

    // Lock the face mutex first
    mesh.entityMutex(2).lock();

    // Signal the calling thread
    thread->sendSignal(meshHandler::START);

    // Recast the pointers for the argument
    faceList& faces =
    (
        *(static_cast<faceList*>(thread->operator()(0)))
    );

    labelList& owner =
    (
        *(static_cast<labelList*>(thread->operator()(1)))
    );

    labelList& neighbour =
    (
        *(static_cast<labelList*>(thread->operator()(2)))
    );

    labelListList& faceEdges =
    (
        *(static_cast<labelListList*>(thread->operator()(3)))
    );

    labelListList& faceZoneFaceMap =
    (
        *(static_cast<labelListList*>(thread->operator()(4)))
    );

    // Reorder the faces
    mesh.reOrderFaces
    (
        faces,
        owner,
        neighbour,
        faceEdges,
        faceZoneFaceMap,
        true
    );
}


// Reorder & renumber cells with bandwidth reduction after a topology change
void dynamicTopoFvMesh::reOrderCells
(
    labelListList& cellZoneMap,
    bool threaded
)
{
    // *** Cell renumbering *** //
    // If cells were deleted during topology change, the numerical order ceases
    // to be continuous. Also, cells are always added at the end of the list by
    // virtue of the append method. Thus, cells would now have to be
    // reordered so that bandwidth is reduced and renumbered to be sequential.

    if (debug)
    {
        if (threaded)
        {
            Info << "Thread: " << self() << ": ";
        }

        Info << "ReOrdering cells..." << flush;
    }

    // Allocate for mapping information
    cellMap_.setSize(nCells_, -1);

    label cellInOrder = 0;

    addedCellRenumbering_.clear();

    // Make a copy of the old cell-based lists, and clear them
    label allCells = cells_.size();

    cellList oldCells(allCells);

    forAll(cells_, cellI)
    {
        oldCells[cellI].transfer(cells_[cellI]);
    }

    cells_.setSize(nCells_);

    if (bandWidthReduction_)
    {
        label currentCell;
        SLList<label> nextCell;
        labelList ncc(allCells, 0), visited(allCells, 0);
        labelListList cellCellAddr(allCells);

        // Build a cell-cell addressing list
        forAll(owner_, faceI)
        {
            if ((neighbour_[faceI] > -1) && (owner_[faceI] > -1))
            {
                ncc[owner_[faceI]]++;
                ncc[neighbour_[faceI]]++;
            }
        }

        forAll(cellCellAddr, cellI)
        {
            cellCellAddr[cellI].setSize(ncc[cellI]);

            // Mark off deleted cells as "visited"
            if (ncc[cellI] == 0)
            {
                visited[cellI] = 1;
            }
        }

        ncc = 0;

        forAll(owner_, faceI)
        {
            if ((owner_[faceI] > -1) && (neighbour_[faceI] > -1))
            {
                cellCellAddr[owner_[faceI]][ncc[owner_[faceI]]++] =
                (
                    neighbour_[faceI]
                );

                cellCellAddr[neighbour_[faceI]][ncc[neighbour_[faceI]]++] =
                (
                    owner_[faceI]
                );
            }
        }

        // Let's get to the "business bit" of the band-compression
        forAll(visited, cellI)
        {
            // Find the first cell that has not been visited yet
            if (visited[cellI] == 0)
            {
                // Use this cell as a start
                currentCell = cellI;

                nextCell.append(currentCell);

                // Loop through the nextCell list. Add the first cell
                // into the cell order if it has not already been visited
                // and ask for its neighbours. If the neighbour in question
                // has not been visited, add it to the end of the nextCell list
                while (nextCell.size() > 0)
                {
                    currentCell = nextCell.removeHead();

                    if (visited[currentCell] == 0)
                    {
                        // Mark as visited and update cell mapping info
                        visited[currentCell] = 1;

                        if (currentCell < nOldCells_)
                        {
                            cellMap_[cellInOrder] = currentCell;
                            reverseCellMap_[currentCell] = cellInOrder;
                        }
                        else
                        {
                            addedCellRenumbering_.insert
                            (
                                currentCell,
                                cellInOrder
                            );
                        }

                        // Insert entities into local lists...
                        cells_[cellInOrder].transfer(oldCells[currentCell]);

                        cellInOrder++;

                        // Find if the neighbours have been visited
                        const labelList& neighbours = cellCellAddr[currentCell];

                        forAll(neighbours, nI)
                        {
                            if (visited[neighbours[nI]] == 0)
                            {
                                // Not visited, add to the list
                                nextCell.append(neighbours[nI]);
                            }
                        }
                    }
                }
            }
        }

        if (debug > 1)
        {
            if (sum(visited) != allCells)
            {
                FatalErrorIn("dynamicTopoFvMesh::reOrderCells()") << nl
                    << " Algorithm did not visit every cell in the mesh."
                    << " Something's messed up." << nl
                    << abort(FatalError);
            }
        }
    }
    else
    {
        // No bandwidth reduction. Fill-in sequentially.
        for (label cellI = 0; cellI < nOldCells_; cellI++)
        {
            // Check if this is a deleted cell
            if (reverseCellMap_[cellI] == -1)
            {
                continue;
            }

            // Update the cell info
            cells_[cellInOrder].transfer(oldCells[cellI]);

            // Update maps
            cellMap_[cellInOrder] = cellI;
            reverseCellMap_[cellI] = cellInOrder;

            // Update the counter
            cellInOrder++;
        }

        for (label cellI = nOldCells_; cellI < allCells; cellI++)
        {
            // Was this cell removed after addition?
            if (deletedCells_.found(cellI))
            {
                continue;
            }

            // Update the cell info
            cells_[cellInOrder].transfer(oldCells[cellI]);

            // Put inserted cells in a seperate hashSet
            addedCellRenumbering_.insert(cellI, cellInOrder);

            // Update the counter
            cellInOrder++;
        }

        // Final check to ensure everything went okay
        if (cellInOrder != nCells_)
        {
            FatalErrorIn("dynamicTopoFvMesh::reOrderCells()") << nl
                << " Algorithm did not visit every cell in the mesh."
                << " Something's messed up." << nl
                << abort(FatalError);
        }
    }

    // Now that we're done preparing the cell maps, unlock the cell mutex
    if (threaded)
    {
        entityMutex_[3].unlock();
    }

    // Renumber all maps.
    forAll(cellsFromPoints_, cellI)
    {
        objectMap& thisMap = cellsFromPoints_[cellI];

        if (thisMap.index() < nOldCells_)
        {
            thisMap.index() = reverseCellMap_[thisMap.index()];
        }
        else
        {
            thisMap.index() = addedCellRenumbering_[thisMap.index()];
        }
    }

    forAll(cellsFromEdges_, cellI)
    {
        objectMap& thisMap = cellsFromEdges_[cellI];

        if (thisMap.index() < nOldCells_)
        {
            thisMap.index() = reverseCellMap_[thisMap.index()];
        }
        else
        {
            thisMap.index() = addedCellRenumbering_[thisMap.index()];
        }
    }

    forAll(cellsFromFaces_, cellI)
    {
        objectMap& thisMap = cellsFromFaces_[cellI];

        if (thisMap.index() < nOldCells_)
        {
            thisMap.index() = reverseCellMap_[thisMap.index()];
        }
        else
        {
            thisMap.index() = addedCellRenumbering_[thisMap.index()];
        }
    }

    forAll(cellsFromCells_, cellI)
    {
        objectMap& thisMap = cellsFromCells_[cellI];

        if (thisMap.index() < nOldCells_)
        {
            thisMap.index() = reverseCellMap_[thisMap.index()];
        }
        else
        {
            thisMap.index() = addedCellRenumbering_[thisMap.index()];
        }
    }

    // Prepare the cellZoneMap
    cellZoneMesh& cellZones = polyMesh::cellZones();

    labelListList newCellZoneAddr(cellZones.size());
    labelList nCellsInZone(cellZones.size(), 0);

    // Prepare zone maps.
    forAll(cellZones, czI)
    {
        // Get the list of old cells
        const labelList& oldAddr = cellZones[czI];

        label& curNCells = nCellsInZone[czI];

        // First count the actual number of cells in each zone.
        forAll(oldAddr, cellI)
        {
            // Was this zone cell deleted? Don't count it.
            if (reverseCellMap_[oldAddr[cellI]] != -1)
            {
                curNCells++;
            }
        }

        // Check for added cells as well
        forAllIter(Map<label>, addedCellZones_, cIter)
        {
            if (cIter() == czI)
            {
                curNCells++;
            }
        }

        label cIndex = 0;
        labelList& newAddr = newCellZoneAddr[czI];

        // Set the sizes first
        newAddr.setSize(curNCells);
        cellZoneMap[czI].setSize(curNCells, -1);

        // Add existing zone cells which have been renumbered.
        forAll(oldAddr, cellI)
        {
            if (reverseCellMap_[oldAddr[cellI]] != -1)
            {
                cellZoneMap[czI][cIndex] = oldAddr[cellI];
                newAddr[cIndex] = reverseCellMap_[oldAddr[cellI]];
                cIndex++;
            }
        }

        // Next, add the newly added zone cells.
        forAllIter(Map<label>, addedCellZones_, cIter)
        {
            if (cIter() == czI)
            {
                newAddr[cIndex++] = addedCellRenumbering_[cIter.key()];
            }
        }

        // Finally, assign addressing to this zone.
        cellZones[czI] = newCellZoneAddr[czI];
    }

    // Reset all zones
    cellZones.updateMesh();

    if (debug)
    {
        Info << "Done." << endl;
    }
}


// Static equivalent for multi-threading
void dynamicTopoFvMesh::reOrderCellsThread
(
    void *argument
)
{
    // Recast the argument
    meshHandler *thread = static_cast<meshHandler*>(argument);

    dynamicTopoFvMesh& mesh = thread->reference();

    // Lock the cell mutex first
    mesh.entityMutex(3).lock();

    // Signal the calling thread
    thread->sendSignal(meshHandler::START);

    // Recast the pointers for the argument
    labelListList& cellZoneMap =
    (
        *(static_cast<labelListList*>(thread->operator()(0)))
    );

    // Reorder the cells
    mesh.reOrderCells(cellZoneMap, true);
}


// Reorder the faces in upper-triangular order, and generate mapping information
void dynamicTopoFvMesh::reOrderMesh
(
    pointField& points,
    pointField& preMotionPoints,
    edgeList& edges,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    labelListList& faceEdges,
    labelListList& edgeFaces,
    labelListList& pointZoneMap,
    labelListList& faceZoneFaceMap,
    labelListList& cellZoneMap
)
{
    if (debug)
    {
        Info << endl;
        Info << "=================" << endl;
        Info << " Mesh reOrdering " << endl;
        Info << "=================" << endl;
        Info << "Mesh Info [n]:" << endl;
        Info << "Points: " << nOldPoints_ << endl;
        Info << "Edges: " << nOldEdges_ << endl;
        Info << "Faces: " << nOldFaces_ << endl;
        Info << "Cells: " << nOldCells_ << endl;
        Info << "Internal Edges: " << nOldInternalEdges_ << endl;
        Info << "Internal Faces: " << nOldInternalFaces_ << endl;
        if (debug > 1)
        {
            Info << "Patch Starts [Face]: " << oldPatchStarts_ << endl;
            Info << "Patch Sizes  [Face]: " << oldPatchSizes_ << endl;
            Info << "Patch Starts [Edge]: " << oldEdgePatchStarts_ << endl;
            Info << "Patch Sizes  [Edge]: " << oldEdgePatchSizes_ << endl;
        }
        Info << "=================" << endl;
        Info << "Mesh Info [n+1]:" << endl;
        Info << "Points: " << nPoints_ << endl;
        Info << "Edges: " << nEdges_ << endl;
        Info << "Faces: " << nFaces_ << endl;
        Info << "Cells: " << nCells_ << endl;
        Info << "Internal Edges: " << nInternalEdges_ << endl;
        Info << "Internal Faces: " << nInternalFaces_ << endl;
        if (debug > 1)
        {
            Info << "Patch Starts [Face]: " << patchStarts_ << endl;
            Info << "Patch Sizes: [Face]: " << patchSizes_ << endl;
            Info << "Patch Starts [Edge]: " << edgePatchStarts_ << endl;
            Info << "Patch Sizes: [Edge]: " << edgePatchSizes_ << endl;
        }
        Info << "=================" << endl;

        // Check connectivity structures for consistency
        // before entering the reOrdering phase.
        checkConnectivity();
    }

    if (threader_->multiThreaded())
    {
        // Initialize multi-threaded reOrdering
        threadedMeshReOrdering
        (
            points,
            preMotionPoints,
            edges,
            faces,
            owner,
            neighbour,
            faceEdges,
            edgeFaces,
            pointZoneMap,
            faceZoneFaceMap,
            cellZoneMap
        );
    }
    else
    {
        // Reorder the points
        reOrderPoints(points, preMotionPoints, pointZoneMap);

        // Reorder the cells
        reOrderCells(cellZoneMap);

        // Reorder the faces
        reOrderFaces(faces, owner, neighbour, faceEdges, faceZoneFaceMap);

        // Reorder the edges
        reOrderEdges(edges, edgeFaces, faceEdges);
    }
}


// Invoke reOrdering with multiple threads
void dynamicTopoFvMesh::threadedMeshReOrdering
(
    pointField& points,
    pointField& preMotionPoints,
    edgeList& edges,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    labelListList& faceEdges,
    labelListList& edgeFaces,
    labelListList& pointZoneMap,
    labelListList& faceZoneFaceMap,
    labelListList& cellZoneMap
)
{
    // For reOrdering, one handler for each reOrdering method
    PtrList<meshHandler> reOrderPtr(4);

    // Initialize reOrdering handlers
    forAll(reOrderPtr, memberI)
    {
        reOrderPtr.set
        (
            memberI,
            new meshHandler(*this, threader())
        );
    }

    // Set argument sizes for individual members

    // Points take three arguments
    // (Two pointFields and one labelListList)
    reOrderPtr[0].setSize(3);

    // Prepare pointers for point reOrdering
    reOrderPtr[0].set(0, &points);
    reOrderPtr[0].set(1, &preMotionPoints);
    reOrderPtr[0].set(2, &pointZoneMap);

    // Edges take three arguments
    // (One edgeList and two labelListLists)
    reOrderPtr[1].setSize(3);

    // Prepare pointers for edge reOrdering
    reOrderPtr[1].set(0, &edges);
    reOrderPtr[1].set(1, &edgeFaces);
    reOrderPtr[1].set(2, &faceEdges);

    // Faces take five arguments
    // (One faceList, two labelLists, and two labelListLists)
    reOrderPtr[2].setSize(5);

    // Prepare pointers for face reOrdering
    reOrderPtr[2].set(0, &faces);
    reOrderPtr[2].set(1, &owner);
    reOrderPtr[2].set(2, &neighbour);
    reOrderPtr[2].set(3, &faceEdges);
    reOrderPtr[2].set(4, &faceZoneFaceMap);

    // Cells take one argument
    // (One labelListList)
    reOrderPtr[3].setSize(1);

    // Prepare pointers for cell reOrdering
    reOrderPtr[3].set(0, &cellZoneMap);

    // Set the thread scheduling sequence
    labelList reOrderSeq(4, -1);

    // Points, cells, faces and edges (in that order)
    reOrderSeq[0] = 0;
    reOrderSeq[1] = 3;
    reOrderSeq[2] = 2;
    reOrderSeq[3] = 1;

    // Lock all slave threads first
    forAll(reOrderSeq, i)
    {
        reOrderPtr[reOrderSeq[i]].lock(meshHandler::START);
        reOrderPtr[reOrderSeq[i]].unsetPredicate(meshHandler::START);
    }

    // Submit points to the work queue
    threader_->addToWorkQueue
    (
        &reOrderPointsThread,
        &(reOrderPtr[0])
    );

    // Wait for a signal from this thread before moving on.
    reOrderPtr[0].waitForSignal(meshHandler::START);

    // Submit cells to the work queue
    threader_->addToWorkQueue
    (
        &reOrderCellsThread,
        &(reOrderPtr[3])
    );

    // Wait for a signal from this thread before moving on.
    reOrderPtr[3].waitForSignal(meshHandler::START);

    // Submit faces to the work queue
    threader_->addToWorkQueue
    (
        &reOrderFacesThread,
        &(reOrderPtr[2])
    );

    // Wait for a signal from this thread before moving on.
    reOrderPtr[2].waitForSignal(meshHandler::START);

    // Lock the edge stop-mutex
    reOrderPtr[1].lock(meshHandler::STOP);
    reOrderPtr[1].unsetPredicate(meshHandler::STOP);

    // Submit edges to the work queue
    threader_->addToWorkQueue
    (
        &reOrderEdgesThread,
        &(reOrderPtr[1])
    );

    // Wait for a signal from this thread before moving on.
    reOrderPtr[1].waitForSignal(meshHandler::START);

    // Wait for edges to be reOrdered before moving on.
    reOrderPtr[1].waitForSignal(meshHandler::STOP);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
