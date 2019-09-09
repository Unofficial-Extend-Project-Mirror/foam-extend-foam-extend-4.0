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

\*---------------------------------------------------------------------------*/

#include "removeFaces.H"
#include "polyMesh.H"
#include "directTopoChange.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "polyModifyFace.H"
#include "polyRemoveFace.H"
#include "polyRemoveCell.H"
#include "polyRemovePoint.H"
#include "syncTools.H"
#include "OFstream.H"
#include "indirectPrimitivePatch.H"
#include "foamTime.H"
#include "faceSet.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class TopoChangeEngine>
void Foam::removeFaces::mergeFaces
(
    const labelList& cellRegion,
    const labelList& cellRegionMaster,
    const labelHashSet& pointsToRemove,
    const labelList& faceLabels,

    TopoChangeEngine& ref
) const
{
    // Construct addressing engine from faceLabels (in order of faceLabels as
    // well)
    indirectPrimitivePatch fp
    (
        IndirectList<face>
        (
            mesh_.faces(),
            faceLabels
        ),
        mesh_.points()
    );

    // Get outside vertices (in local vertex numbering)
    if (fp.edgeLoops().size() != 1)
    {
        writeOBJ(fp, mesh_.time().path()/"facesToBeMerged.obj");

        FatalErrorInFunction
            << "Cannot merge faces " << faceLabels
            << " into single face since outside vertices " << fp.edgeLoops()
            << " do not form single loop but form " << fp.edgeLoops().size()
            << " loops instead." << abort(FatalError);
    }

    const labelList& edgeLoop = fp.edgeLoops()[0];

    // Get outside vertices in order of one of the faces in faceLabels.
    // (this becomes the master face)
    // Find the first face that uses edgeLoop[0] and edgeLoop[1] as consecutive
    // vertices.

    label masterIndex = -1;
    bool reverseLoop = false;

    const labelList& pFaces = fp.pointFaces()[edgeLoop[0]];

    // Find face among pFaces which uses edgeLoop[1]
    forAll(pFaces, i)
    {
        const label& faceI = pFaces[i];

        const face& f = fp.localFaces()[faceI];

        label index1 = findIndex(f, edgeLoop[1]);

        if (index1 != -1)
        {
            // Check whether consecutive to edgeLoop[0]
            label index0 = findIndex(f, edgeLoop[0]);

            if (index0 != -1)
            {
                if (index1 == f.fcIndex(index0))
                {
                    masterIndex = faceI;
                    reverseLoop = false;
                    break;
                }
                else if (index1 == f.rcIndex(index0))
                {
                    masterIndex = faceI;
                    reverseLoop = true;
                    break;
                }
            }
        }
    }

    if (masterIndex == -1)
    {
        writeOBJ(fp, mesh_.time().path()/"facesToBeMerged.obj");

        FatalErrorInFunction
            << "Problem." << abort(FatalError);
    }


    // Modify the master face

    // Modify first face
    label faceI = faceLabels[masterIndex];

    label own = mesh_.faceOwner()[faceI];

    if (cellRegion[own] != -1)
    {
        own = cellRegionMaster[cellRegion[own]];
    }

    // Set face information
    label patchID, zoneID, zoneFlip;
    meshTools::setFaceInfo(mesh_, faceI, patchID, zoneID, zoneFlip);

    label nei = -1;

    if (mesh_.isInternalFace(faceI))
    {
        nei = mesh_.faceNeighbour()[faceI];

        if (cellRegion[nei] != -1)
        {
            nei = cellRegionMaster[cellRegion[nei]];
        }
    }

    // Collect non-removed face vertices
    dynamicLabelList faceVerts(edgeLoop.size());

    // Get mesh points
    const labelList& meshPoints = fp.meshPoints();

    forAll(edgeLoop, i)
    {
        const label& pointI = meshPoints[edgeLoop[i]];

        if (pointsToRemove.found(pointI))
        {
            // Point should be removed, nothing to do
        }
        else
        {
            faceVerts.append(pointI);
        }
    }

    face mergedFace;
    mergedFace.transfer(faceVerts);
    faceVerts.clear();

    if (reverseLoop)
    {
        // Reverse the merged face
        reverse(mergedFace);
    }

    // Finally modify merged face
    modifyFace
    (
        mergedFace,         // modified face
        faceI,              // label of face being modified
        own,                // owner
        nei,                // neighbour
        false,              // face flip
        patchID,            // patch for face
        false,              // remove from zone
        zoneID,             // zone for face
        zoneFlip,           // face flip in zone

        ref                 // topo change engine
    );


    // Remove all but master face
    forAll(faceLabels, patchFaceI)
    {
        if (patchFaceI != masterIndex)
        {
            ref.setAction(polyRemoveFace(faceLabels[patchFaceI], faceI));
        }
    }
}


template<class TopoChangeEngine>
void Foam::removeFaces::modifyFace
(
    const face& f,
    const label masterFaceID,
    const label own,
    const label nei,
    const bool flipFaceFlux,
    const label newPatchID,
    const bool removeFromZone,
    const label zoneID,
    const bool zoneFlip,

    TopoChangeEngine& ref
) const
{
    // Print info with deep debug level only
    if (debug > 1)
    {
        if (mesh_.isInternalFace(masterFaceID))
        {
            Pout<< "Modifying face " << masterFaceID
                << ", old verts: " << mesh_.faces()[masterFaceID]
                << " for new verts:" << f
                << nl
                << " or for new owner " << own
                << " (old owner: " << mesh_.faceOwner()[masterFaceID] << ")"
                << nl
                << " or for new nei " << nei << " (old neighbour: "
                << mesh_.faceNeighbour()[masterFaceID] << ")"
                << endl;
        }
    }

    if ((nei == -1) || (own < nei))
    {
        // No need to revert the face, use original owner/neighbour
        ref.setAction
        (
            polyModifyFace
            (
                f,              // modified face
                masterFaceID,   // label of face being modified
                own,            // owner
                nei,            // neighbour
                flipFaceFlux,   // face flip
                newPatchID,     // patch for face
                removeFromZone, // remove from zone
                zoneID,         // zone for face
                zoneFlip        // face flip in zone
            )
        );
    }
    else
    {
        // Revert the face and flip owner/neighbour
        ref.setAction
        (
            polyModifyFace
            (
                f.reverseFace(),// modified face
                masterFaceID,   // label of face being modified
                nei,            // owner
                own,            // neighbour
                flipFaceFlux,   // face flip
                newPatchID,     // patch for face
                removeFromZone, // remove from zone
                zoneID,         // zone for face
                zoneFlip        // face flip in zone
            )
        );
    }
}


template<class TopoChangeEngine>
void Foam::removeFaces::setRefinement
(
    const labelList& faceLabels,
    const labelList& cellRegion,
    const labelList& pointRegionMaster,
    labelList& cellRegionMaster,

    TopoChangeEngine& ref
) const
{
    if (debug)
    {
        const faceSet facesToRemove
        (
            mesh_,
            "facesToRemove",
            labelHashSet(faceLabels)
        );

        Pout<< "Writing faces to remove to faceSet " << facesToRemove.name()
            << endl;

        facesToRemove.write();
    }

    // Get number of mesh faces
    const label nFaces = mesh_.nFaces();

    // Mark-up field for all faces that need to be removed
    boolList removedFace(nFaces, false);
    forAll(faceLabels, i)
    {
        const label& faceI = faceLabels[i];

        if (!mesh_.isInternalFace(faceI))
        {
            FatalErrorInFunction
                << "Face " << faceI << " is not an internal faces, therefore"
                << " it cannot be removed. Check faceLabels argument."
                << abort(FatalError);
        }

        removedFace[faceI] = true;
    }


    // PART 1: Collect edges to be removed
    labelHashSet edgesToRemove(faceLabels.size());

    // Region for each face:
    // -1 = removed faces
    // -2 = regions consisting of single face only
    labelList faceRegion(nFaces, -1);

    // Number of connected face regions
    label nRegions = 0;

    {
        // Number of edges per non-removed face. See below regarding
        // initialization
        labelList nFacesPerEdge(mesh_.nEdges(), -1);

        // Get necessary mesh data
        const labelListList& meshFaceEdges = mesh_.faceEdges();
        const labelListList& meshEdgeFaces = mesh_.edgeFaces();
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        // Count usage of edges by non-removed faces.
        forAll(faceLabels, i)
        {
            // Get face Index and edges of this face
            const label& faceI = faceLabels[i];
            const labelList& fEdges = meshFaceEdges[faceI];

            forAll(fEdges, j)
            {
                // Get edge index
                const label& edgeI = fEdges[j];

                if (nFacesPerEdge[edgeI] == -1)
                {
                    // Number of faces for this edge is not set, set to size - 1
                    nFacesPerEdge[edgeI] = meshEdgeFaces[edgeI].size() - 1;
                }
                else
                {
                    // Decrement number of faces for this edge
                    --nFacesPerEdge[edgeI];
                }
            }
        }

        // Count usage for edges not on faces-to-be-removed.
        // Note that this only needs to be done for possibly coupled edges
        // so we could choose to loop only over boundary faces and use faceEdges
        // of those.
        forAll(meshEdgeFaces, edgeI)
        {
            if (nFacesPerEdge[edgeI] == -1)
            {
                // Edge not yet handled in loop above so is not used by any
                // face to be removed. Get edge faces
                const labelList& eFaces = meshEdgeFaces[edgeI];

                // Number of faces for this edge is greater than 2, set size
                if (eFaces.size() > 2)
                {
                    nFacesPerEdge[edgeI] = eFaces.size();
                }
                else if (eFaces.size() == 2)
                {
                    // nFacesPerEdge already -1 so do nothing
                }
                else
                {
                    // Error: edge has less than two faces. Print additional
                    // information and issue an error

                    // Get edge and all its cells
                    const edge& e = mesh_.edges()[edgeI];
                    const labelListList& ec = mesh_.edgeCells();

                    // Get point cells
                    const labelListList& pc = mesh_.pointCells();

                    // Write mesh before termination
                    mesh_.write();

                    // Write data for debugging
                    FatalErrorInFunction
                        << "Problem : edge has too few face neighbours:"
                        << eFaces << endl
                        << "edge:" << edgeI
                        << " vertices:" << e
                        << " coords:" << mesh_.points()[e[0]]
                        << mesh_.points()[e[1]]
                        << endl
                        << "Edge cells: " << ec[edgeI] << nl
                        << "First point cells: " << pc[e[0]] << nl
                        << "Second point cells: " << pc[e[1]] << nl
                        << abort(FatalError);
                }
            }
        }


        if (debug)
        {
            // Write edges with two faces
            OFstream str(mesh_.time().path()/"edgesWithTwoFaces.obj");

            Pout<< "Writing edgesWithTwoFaces to " << str.name() << endl;

            // Vertex counter
            label vertI = 0;

            // Get mesh points and edges
            const pointField& meshPoints = mesh_.points();
            const edgeList& meshEdges = mesh_.edges();

            forAll(nFacesPerEdge, edgeI)
            {
                if (nFacesPerEdge[edgeI] == 2)
                {
                    // Edge will get removed, write data
                    const edge& e = meshEdges[edgeI];

                    meshTools::writeOBJ(str, meshPoints[e[0]]);
                    ++vertI;
                    meshTools::writeOBJ(str, meshPoints[e[1]]);
                    ++vertI;
                    str<< "l " << vertI - 1 << ' ' << vertI << nl;
                }
            }
        }


        // At this point, all affected edges have the number of unremoved faces

        // Filter for edges inbetween two remaining boundary faces that
        // make too big an angle.
        forAll(nFacesPerEdge, edgeI)
        {
            if (nFacesPerEdge[edgeI] == 2)
            {
                // See if these are two boundary faces
                label f0 = -1;
                label f1 = -1;

                const labelList& eFaces = meshEdgeFaces[edgeI];

                forAll(eFaces, i)
                {
                    const label& faceI = eFaces[i];

                    if (!removedFace[faceI] && !mesh_.isInternalFace(faceI))
                    {
                        if (f0 == -1)
                        {
                            f0 = faceI;
                        }
                        else
                        {
                            f1 = faceI;
                            break;
                        }
                    }
                }

                if (f0 != -1 && f1 != -1)
                {
                    // Edge has two boundary faces remaining, see if they should
                    // be merged.
                    const label patch0 = patches.whichPatch(f0);
                    const label patch1 = patches.whichPatch(f1);

                    if (patch0 != patch1)
                    {
                        // Different patches. Do not merge edge.
                        WarningIn("removeFaces::setRefinement")
                            << "Not merging faces " << f0 << " and "
                            << f1 << " across patch boundary edge " << edgeI
                            << " since they are not on the same patch."
                            << endl;

                        // Set number of faces to 3 to preserve the face
                        nFacesPerEdge[edgeI] = 3;
                    }
                    else if (minCos_ < 1 && minCos_ > -1)
                    {
                        const polyPatch& pp0 = patches[patch0];
                        const vectorField& n0 = pp0.faceNormals();

                        if
                        (
                            mag
                            (
                                n0[f0 - pp0.start()]
                              & n0[f1 - pp0.start()]
                            )
                            < minCos_
                        )
                        {
                            WarningIn("removeFaces::setRefinement")
                                << "Not merging faces " << f0 << " and "
                                << f1 << " across edge " << edgeI
                                << " since the angle between them is too large."
                                << endl;

                            // Set number of faces to 3 to preserve the face
                            nFacesPerEdge[edgeI] = 3;
                        }
                    }
                }
                else if (f0 != -1 || f1 != -1)
                {
                    const edge& e = mesh_.edges()[edgeI];

                    // Only found one boundary face. Problem
                    FatalErrorInFunction
                        << "Problem : edge would have one boundary face"
                        << " and one internal face using it." << endl
                        << "Your remove pattern is probably incorrect." << endl
                        << "edge:" << edgeI
                        << " nFaces:" << nFacesPerEdge[edgeI]
                        << " vertices:" << e
                        << " coords:" << mesh_.points()[e[0]]
                        << mesh_.points()[e[1]]
                        << " face0:" << f0
                        << " face1:" << f1
                        << abort(FatalError);
                }
            }
        }


        // Check locally (before synchronizing) for strangeness
        forAll(nFacesPerEdge, edgeI)
        {
            if (nFacesPerEdge[edgeI] == 1)
            {
                const edge& e = mesh_.edges()[edgeI];

                FatalErrorInFunction
                    << "Problem : edge would get 1 face using it only"
                    << " edge:" << edgeI
                    << " nFaces:" << nFacesPerEdge[edgeI]
                    << " vertices:" << e
                    << " coords:" << mesh_.points()[e[0]]
                    << ' ' << mesh_.points()[e[1]]
                    << abort(FatalError);
            }

            // Could check here for boundary edge with <= 1 faces remaining
        }


        // Synchronize edge usage. This is to make sure that both sides remove
        // (or not remove) an edge on the boundary at the same time.
        //
        // Coupled edges (edge0, edge1 are opposite each other)
        // a. edge not on face to be removed, edge has >= 3 faces
        // b.  ,,                             edge has 2 faces
        // c. edge has >= 3 remaining faces
        // d. edge has 2 remaining faces (assume angle>minCos already handled)
        //
        // - a + a: do not remove edge
        // - a + b: do not remove edge
        // - a + c: do not remove edge
        // - a + d: do not remove edge
        //
        // - b + b: do not remove edge
        // - b + c: do not remove edge
        // - b + d: remove edge
        //
        // - c + c: do not remove edge
        // - c + d: do not remove edge
        // - d + d: remove edge
        //
        //
        // So code situation a. with >= 3
        //                   b. with -1
        //                   c. with >=3
        //                   d. with 2
        // then do max and check result.
        //
        // a+a : max(3,3) = 3. do not remove
        // a+b : max(3,-1) = 3. do not remove
        // a+c : max(3,3) = 3. do not remove
        // a+d : max(3,2) = 3. do not remove
        //
        // b+b : max(-1,-1) = -1. do not remove
        // b+c : max(-1,3) = 3. do not remove
        // b+d : max(-1,2) = 2. remove
        //
        // c+c : max(3,3) = 3. do not remove
        // c+d : max(3,2) = 3. do not remove
        //
        // d+d : max(2,2) = 2. remove

        syncTools::syncEdgeList
        (
            mesh_,
            nFacesPerEdge,
            maxEqOp<label>(),
            labelMin,               // guaranteed to be overridden by maxEqOp
            false                   // no separation
        );

        // Convert to labelHashSet
        forAll(nFacesPerEdge, edgeI)
        {
            if (nFacesPerEdge[edgeI] == 0)
            {
                // 0: edge not used anymore.
                edgesToRemove.insert(edgeI);
            }
            else if (nFacesPerEdge[edgeI] == 1)
            {
                // 1: illegal. Tested above.
            }
            else if (nFacesPerEdge[edgeI] == 2)
            {
                // 2: merge faces.
                edgesToRemove.insert(edgeI);
            }
        }

        if (debug)
        {
            // Write edges to remove
            OFstream str(mesh_.time().path()/"edgesToRemove.obj");

            Pout<< "Writing edgesToRemove to " << str.name() << endl;

            // Vertex counter
            label vertI = 0;

            // Get mesh points
            const pointField& meshPoints = mesh_.points();

            forAllConstIter(labelHashSet, edgesToRemove, iter)
            {
                // Edge will get removed.
                const edge& e = mesh_.edges()[iter.key()];

                meshTools::writeOBJ(str, meshPoints[e[0]]);
                ++vertI;
                meshTools::writeOBJ(str, meshPoints[e[1]]);
                ++vertI;
                str<< "l " << vertI - 1 << ' ' << vertI << nl;
            }
        }

        // Walk to fill faceRegion with faces that will be connected across
        // edges that will be removed

        label startFaceI = 0;

        while (true)
        {
            // Find unset region
            for (; startFaceI < nFaces; ++startFaceI)
            {
                if (faceRegion[startFaceI] == -1 && !removedFace[startFaceI])
                {
                    break;
                }
            }

            if (startFaceI == nFaces)
            {
                break;
            }

            // Start walking face-edge-face, crossing edges that will get
            // removed. Every thus connected region will get single region
            // number.
            label nRegion = changeFaceRegion
            (
                cellRegion,
                removedFace,
                nFacesPerEdge,
                startFaceI,
                nRegions,

                faceRegion
            );

            if (nRegion < 1)
            {
                FatalErrorInFunction
                    << "Problem with region number." << abort(FatalError);
            }
            else if (nRegion == 1)
            {
                // Reset face to be single region
                faceRegion[startFaceI] = -2;
            }
            else
            {
                ++nRegions;
            }
        }

        // Check we're deciding the same on both sides. Since the regioning
        // is done based on nFacesPerEdge (which is synced) this should
        // indeed be the case.

        // Create a copy of face region list and swap it
        labelList nbrFaceRegion(faceRegion);
        syncTools::swapFaceList
        (
            mesh_,
            nbrFaceRegion,
            false                   // no separation
        );

        // Store data from neighbouring region
        labelList toNbrRegion(nRegions, -1);

        // Get number of internal faces
        const label nInternalFaces = mesh_.nInternalFaces();

        for
        (
            label faceI = nInternalFaces;
            faceI < nFaces;
            ++faceI
        )
        {
            // Get the neighbouring region
            const label& nbrRegion = nbrFaceRegion[faceI];
            const label& myRegion = faceRegion[faceI];

            if (myRegion <= -1 || nbrRegion <= -1)
            {
                if (nbrRegion != myRegion)
                {
                    FatalErrorInFunction
                        << "Inconsistent face region across coupled patches."
                        << endl
                        << "This side has for faceI:" << faceI
                        << " region:" << myRegion << endl
                        << "The other side has region:" << nbrRegion
                        << endl
                        << "(region - 1 means face is to be deleted)"
                        << abort(FatalError);
                }
            }
            else if (toNbrRegion[myRegion] == -1)
            {
                // First visit of region. Store correspondence.
                toNbrRegion[myRegion] = nbrRegion;
            }
            else
            {
                // Second visit of this region
                if (toNbrRegion[myRegion] != nbrRegion)
                {
                    FatalErrorInFunction
                        << "Inconsistent face region across coupled patches."
                        << endl
                        << "This side has for faceI:" << faceI
                        << " region:" << myRegion
                        << " with coupled neighbouring regions:"
                        << toNbrRegion[myRegion] << " and "
                        << nbrRegion
                        << abort(FatalError);
                }
            }
        }
    } // End memory management


    // PART 2: Collect points to remove

    // Create a hash set of points to remove, assuming that each face has 4
    // points to prevent excessive resizing
    labelHashSet pointsToRemove(4*faceLabels.size());

    // Memory management
    {
        // For each point, count the number of edges that will be kept
        // (unremoved edges). Store the ones that are only used by exactly two
        // unremoved edges.

        // Get necessary mesh data
        const label nPoints = mesh_.nPoints();
        const labelListList& meshPointEdges = mesh_.pointEdges();
        const edgeList& meshEdges = mesh_.edges();

        // List containing number of nonremoved edges for each point
        labelList nEdgesPerPoint(nPoints);

        // Initialise to number of edges per point
        forAll(meshPointEdges, pointI)
        {
            nEdgesPerPoint[pointI] = meshPointEdges[pointI].size();
        }

        // Loop through edges to remove
        forAllConstIter(labelHashSet, edgesToRemove, iter)
        {
            // Get edge to be removed
            const edge& e = meshEdges[iter.key()];

            // Loop through both points of the edge and decrement the counter of
            // unremoved edges per point
            forAll(e, i)
            {
                --nEdgesPerPoint[e[i]];
            }
        }

        // Check locally (before synchronizing) for strangeness
        forAll(nEdgesPerPoint, pointI)
        {
            if (nEdgesPerPoint[pointI] == 1)
            {
                FatalErrorInFunction
                    << "Problem : point would get 1 edge using it only."
                    << " pointI:" << pointI
                    << " coord:" << mesh_.points()[pointI]
                    << abort(FatalError);
            }
        }

        // Synchronize point usage. This is to make sure that both sides remove
        // (or don't remove) a point at the boundary at the same time.
        syncTools::syncPointList
        (
            mesh_,
            nEdgesPerPoint,
            maxEqOp<label>(),
            labelMin,
            false                   // no separation
        );

        forAll(nEdgesPerPoint, pointI)
        {
            if (nEdgesPerPoint[pointI] == 0)
            {
                pointsToRemove.insert(pointI);
            }
            else if (nEdgesPerPoint[pointI] == 1)
            {
                // Already checked before
            }
            else if (nEdgesPerPoint[pointI] == 2)
            {
                // Remove point and merge edges
                pointsToRemove.insert(pointI);
            }
        }
    }


    if (debug)
    {
        // Write points to remove
        OFstream str(mesh_.time().path()/"pointsToRemove.obj");
        Pout<< "Writing pointsToRemove to " << str.name() << endl;

        // Get mesh points
        const pointField& meshPoints = mesh_.points();

        forAllConstIter(labelHashSet, pointsToRemove, iter)
        {
            meshTools::writeOBJ(str, meshPoints[iter.key()]);
        }
    }


    // PART 3: Collect all faces affected in any way by removal of points,
    // edges, faces and cells

    // Note: return type of affectedFaces is Xfer<boolList>, so there is no
    // unnecessary copying
    boolList affectedFace
    (
        affectedFaces
        (
            cellRegion,
            cellRegionMaster,
            faceLabels,
            edgesToRemove,
            pointsToRemove
        )
    );

    // Now the data is complete:
    // - faceLabels         : faces to remove (sync since no boundary faces)
    // - cellRegion/Master  : cells to remove (sync since cells)
    // - pointsToRemove     : points to remove (sync)
    // - faceRegion         : connected face region of faces to be merged (sync)
    // - affectedFace       : faces with points removed and/or owner/neighbour
    //                        changed (non sync)
    // We can start inserting mesh modifier instructions and keep track of
    // changed faces.


    // PART 4: Do all removals first

    // Remove split faces
    forAll(faceLabels, labelI)
    {
        // Get face index
        const label& faceI = faceLabels[labelI];

        // Remove face if not yet uptodate (which can't happen here; but want to
        // be consistent with rest of face removals/modifications)
        if (affectedFace[faceI])
        {
            // Mark face as unaffected
            affectedFace[faceI] = false;

            // Insert face removal instruction into topo change engine
            ref.setAction(polyRemoveFace(faceI, -1));
        }
    }

    // Remove points
    forAllConstIter(labelHashSet, pointsToRemove, iter)
    {
        // Get point index
        const label& pointI = iter.key();

        // Inser point removal instruction into topo change engine
        ref.setAction(polyRemovePoint(pointI, -1));
    }

    // Add master cells for correct mapping
    forAll (cellRegionMaster, regionI)
    {
        // Note: it is legal to have cellRegionMaster = -1 if the region
        // has been created and them marged into another region.
        // Such masters will also have pointRegionMaster = -1 and should
        // be ignored.  HJ, 6/Sep/2019

        // Additionally protect for old directTopoChangers which do not
        // identify points for mapping.  Non-existent pointRegionMaster
        // is rejected
        if (cellRegionMaster[regionI] > -1 && pointRegionMaster[regionI] > -1)
        {
            // Add master cell from master point for correct mapping
            cellRegionMaster[regionI] =
                ref.setAction
                (
                    polyAddCell
                    (
                        pointRegionMaster[regionI], // masterPointID
                        -1,                         // masterEdgeID
                        -1,                         // masterFaceID
                        -1,                         // masterCellID
                        mesh_.cellZones().whichZone(cellRegionMaster[regionI])
                    )
                );
        }
    }

    // Remove cells
    forAll(cellRegion, cellI)
    {
        label region = cellRegion[cellI];

        // Old check is acceptable: for mapping from point, the cellRegionMaster
        // has been replaced in polyAddPoint
        // HJ, 6/Sep/2019
        if (region != -1 && (cellI != cellRegionMaster[region]))
        {
            ref.setAction(polyRemoveCell(cellI, cellRegionMaster[region]));
        }
    }


    // PART 5: Merge faces across edges to be merged

    // Memory management
    {
        // Invert faceRegion so we get region to faces
        labelListList regionToFaces(invertOneToMany(nRegions, faceRegion));

        // Loop through regions
        forAll(regionToFaces, regionI)
        {
            // Get region faces
            const labelList& rFaces = regionToFaces[regionI];

            if (rFaces.size() <= 1)
            {
                FatalErrorInFunction
                    << "Region: " << regionI
                    << " contains only these faces: " << rFaces
                    << abort(FatalError);
            }

            // rFaces[0] is master, rest get removed
            mergeFaces
            (
                cellRegion,
                cellRegionMaster,
                pointsToRemove,
                rFaces,
                ref
            );

            // Mark region faces as visited (non affected anymore)
            forAll(rFaces, i)
            {
                affectedFace[rFaces[i]] = false;
            }
        }
    }


    // PART 6: Remaining affected faces

    // Get necessary mesh data
    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();


    // Check any remaining faces that have not been updated for either new
    // owner/neighbour or points removed
    forAll(affectedFace, faceI)
    {
        if (affectedFace[faceI])
        {
            // Mark face as unaffected
            affectedFace[faceI] = false;

            // Get filtered face: without points that will be removed
            const face f(filterFace(pointsToRemove, faceI));

            // Get owner of the face and change it
            label own = owner[faceI];
            if (cellRegion[own] != -1)
            {
                own = cellRegionMaster[cellRegion[own]];
            }

            // Set face information
            label patchID, zoneID, zoneFlip;
            meshTools::setFaceInfo(mesh_, faceI, patchID, zoneID, zoneFlip);

            // Get neighbour of the face and change it
            label nei = -1;
            if (mesh_.isInternalFace(faceI))
            {
                nei = neighbour[faceI];
                if (cellRegion[nei] != -1)
                {
                    nei = cellRegionMaster[cellRegion[nei]];
                }
            }

            // Finally modify the face
            modifyFace
            (
                f,                  // modified face
                faceI,              // label of face being modified
                own,                // owner
                nei,                // neighbour
                false,              // face flip
                patchID,            // patch for face
                false,              // remove from zone
                zoneID,             // zone for face
                zoneFlip,           // face flip in zone

                ref                 // topo change engine
            );
        }
    }
}


// ************************************************************************* //
