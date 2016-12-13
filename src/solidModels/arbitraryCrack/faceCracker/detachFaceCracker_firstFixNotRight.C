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

#include "faceCracker.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "polyTopoChange.H"
#include "polyTopoChanger.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceCracker::detachFaceCracker
(
    polyTopoChange& ref
) const
{
    if (debug)
    {
        Pout<< "void faceCracker::detachFaceCracker("
            << "polyTopoChange& ref) const "
            << " for object " << name() << " : "
            << "Detaching faces" << endl;
    }

    const polyMesh& mesh = topoChanger().mesh();
    const faceZoneMesh& zoneMesh = mesh.faceZones();

    const primitiveFacePatch& masterFaceLayer =
        zoneMesh[crackZoneID_.index()]();
    const pointField& points = mesh.points();
    const labelListList& meshEdgeFaces = mesh.edgeFaces();

    const labelList& mp = masterFaceLayer.meshPoints();
    const edgeList& zoneLocalEdges = masterFaceLayer.edges();

    const labelList& meshEdges = zoneMesh[crackZoneID_.index()].meshEdges();

    // Create the points

    labelList addedPoints(mp.size(), -1);

    // Go through boundary edges of the master patch.  If all the faces from
    // this patch are internal, mark the points in the addedPoints lookup
    // with their original labels to stop duplication
    label nIntEdges = masterFaceLayer.nInternalEdges();

    for
    (
        label curEdgeID = nIntEdges;
        curEdgeID < meshEdges.size();
        curEdgeID++
    )
    {
        const labelList& curFaces = meshEdgeFaces[meshEdges[curEdgeID]];

        bool edgeIsInternal = true;

        forAll (curFaces, faceI)
        {
            if (!mesh.isInternalFace(curFaces[faceI]))
            {
                // The edge belongs to a boundary face
                edgeIsInternal = false;
                break;
            }
        }

        if (edgeIsInternal)
        {
            // Reset the point creation
            addedPoints[zoneLocalEdges[curEdgeID].start()] =
                mp[zoneLocalEdges[curEdgeID].start()];

            addedPoints[zoneLocalEdges[curEdgeID].end()] =
                mp[zoneLocalEdges[curEdgeID].end()];
        }
    }

    // Create new points for face zone
    forAll (addedPoints, pointI)
    {
        if (addedPoints[pointI] < 0)
        {
            addedPoints[pointI] =
                ref.setAction
                (
                    polyAddPoint
                    (
                        points[mp[pointI]],        // point
                        mp[pointI],                // master point
                        -1,                        // zone ID
                        true                       // supports a cell
                    )
                );
        }
    }

    // Modify faces in the master zone and duplicate for the slave zone

    const labelList& mf = zoneMesh[crackZoneID_.index()];
    const boolList& mfFlip = zoneMesh[crackZoneID_.index()].flipMap();
    const faceList& zoneFaces = masterFaceLayer.localFaces();

    const faceList& faces = mesh.faces();
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    forAll (mf, faceI)
    {
        const label curFaceID = mf[faceI];

        // Build the face for the slave patch by renumbering
        const face oldFace = zoneFaces[faceI].reverseFace();

        face newFace(oldFace.size());

        forAll (oldFace, pointI)
        {
            newFace[pointI] = addedPoints[oldFace[pointI]];
        }

        if (mfFlip[faceI])
        {
            // Face needs to be flipped for the master patch
            ref.setAction
            (
                polyModifyFace
                (
                    faces[curFaceID].reverseFace(), // modified face
                    curFaceID,                  // label of face being modified
                    nei[curFaceID],                 // owner
                    -1,                             // neighbour
                    true,                           // face flip
                    crackPatchID_.index(),          // patch for face
                    false,                          // remove from zone
                    crackZoneID_.index(),           // zone for face
                    !mfFlip[faceI]                  // face flip in zone
                )
            );

            // Add renumbered face into the slave patch
            ref.setAction
            (
                polyAddFace
                (
                    newFace,                        // face
                    own[curFaceID],                 // owner
                    -1,                             // neighbour
                    -1,                             // master point
                    -1,                             // master edge
                    curFaceID,                      // master face
                    false,                          // flip flux
                    crackPatchID_.index(),          // patch to add the face to
                    -1,                             // zone for face
                    false                           // zone flip
                )
            );
        }
        else
        {
            // No flip
            ref.setAction
            (
                polyModifyFace
                (
                    faces[curFaceID],         // modified face
                    curFaceID,                // label of face being modified
                    own[curFaceID],           // owner
                    -1,                       // neighbour
                    false,                    // face flip
                    crackPatchID_.index(),    // patch for face
                    false,                    // remove from zone
                    crackZoneID_.index(),     // zone for face
                    mfFlip[faceI]             // face flip in zone
                )
            );

            // Add renumbered face into the slave patch
            ref.setAction
            (
                polyAddFace
                (
                    newFace,                        // face
                    nei[curFaceID],                 // owner
                    -1,                             // neighbour
                    -1,                             // master point
                    -1,                             // master edge
                    curFaceID,                      // master face
                    true,                           // flip flux
                    crackPatchID_.index(),          // patch to add the face to
                    -1,                             // zone for face
                    false                           // face flip in zone
                )
            );
        }
    }

    // Modify the remaining faces of the master cells to reconnect to the new
    // layer of faces.

    // Algorithm: Go through all the cells of the master zone and make
    // a map of faces to avoid duplicates.  Do not insert the faces in
    // the master patch (as they have already been dealt with).  Make
    // a master layer point renumbering map, which for every point in
    // the master layer gives its new label. Loop through all faces in
    // the map and attempt to renumber them using the master layer
    // point renumbering map.  Once the face is renumbered, compare it
    // with the original face; if they are the same, the face has not
    // changed; if not, modify the face but keep all of its old
    // attributes (apart from the vertex numbers).


    // Create the map of faces in the master cell layer
    // Bug-fix: PC, HJ and ZT, 19 Feb 2013

    labelHashSet masterCellFaceMap(12*mf.size());

    const labelListList& pointFaces = mesh.pointFaces();

    forAll(mf, faceI)
    {
        const labelList& curFacePoints = faces[mf[faceI]];

        forAll(curFacePoints, pointI)
        {
            const labelList& curPointFaces =
                pointFaces[curFacePoints[pointI]];

            forAll(curPointFaces, fI)
            {
                if
                (
                    zoneMesh.whichZone(curPointFaces[fI])
                 != crackZoneID_.index()
                )
                {
                    masterCellFaceMap.insert(curPointFaces[fI]);
                }
            }
        }
    }


//     // Create the map of faces in the master cell layer
//     const labelList& mc =
//         mesh.faceZones()[crackZoneID_.index()].masterCells();

//     labelHashSet masterCellFaceMap(6*mc.size());

//     const cellList& cells = mesh.cells();

//     forAll (mc, cellI)
//     {
//         const labelList& curFaces = cells[mc[cellI]];

//         forAll (curFaces, faceI)
//         {
//             // Check if the face belongs to the master patch; if not add it
//             if (zoneMesh.whichZone(curFaces[faceI]) != crackZoneID_.index())
//             {
//                 masterCellFaceMap.insert(curFaces[faceI]);
//             }
//         }
//     }

//     // Extend the map to include first neighbours of the master cells to
//     // deal with multiple corners.
//     { // Protection and memory management
//         // Make a map of master cells for quick reject
//         labelHashSet mcMap(2*mc.size());

//         forAll (mc, mcI)
//         {
//             mcMap.insert(mc[mcI]);
//         }

//         // Go through all the faces in the masterCellFaceMap.  If the
//         // cells around them are not already used, add all of their
//         // faces to the map
//         const labelList mcf = masterCellFaceMap.toc();

//         forAll (mcf, mcfI)
//         {
//             // Do the owner side
//             const label ownCell = own[mcf[mcfI]];

//             if (!mcMap.found(ownCell))
//             {
//                 // Cell not found. Add its faces to the map
//                 const cell& curFaces = cells[ownCell];

//                 forAll (curFaces, faceI)
//                 {
//                     masterCellFaceMap.insert(curFaces[faceI]);
//                 }
//             }

//             // Do the neighbour side if face is internal
//             if (mesh.isInternalFace(mcf[mcfI]))
//             {
//                 const label neiCell = nei[mcf[mcfI]];

//                 if (!mcMap.found(neiCell))
//                 {
//                     // Cell not found. Add its faces to the map
//                     const cell& curFaces = cells[neiCell];

//                     forAll (curFaces, faceI)
//                     {
//                         masterCellFaceMap.insert(curFaces[faceI]);
//                     }
//                 }
//             }
//         }
//     }

    // Create the master layer point map
    Map<label> masterLayerPointMap(2*mp.size());

    forAll (mp, pointI)
    {
        masterLayerPointMap.insert
        (
            mp[pointI],
            addedPoints[pointI]
        );
    }

    // Grab the list of faces of the master layer
    const labelList masterCellFaces = masterCellFaceMap.toc();

    forAll (masterCellFaces, faceI)
    {
        // Attempt to renumber the face using the masterLayerPointMap.
        // Missing point remain the same

        const label curFaceID = masterCellFaces[faceI];

        const face& oldFace = faces[curFaceID];

        face newFace(oldFace.size());

        bool changed = false;

        forAll (oldFace, pointI)
        {
            if (masterLayerPointMap.found(oldFace[pointI]))
            {
                changed = true;

                newFace[pointI] = masterLayerPointMap.find(oldFace[pointI])();
            }
            else
            {
                newFace[pointI] = oldFace[pointI];
            }
        }

        // If the face has changed, create a modification entry
        if (changed)
        {
            if (mesh.isInternalFace(curFaceID))
            {
                ref.setAction
                (
                    polyModifyFace
                    (
                        newFace,                    // face
                        curFaceID,                  // master face
                        own[curFaceID],             // owner
                        nei[curFaceID],             // neighbour
                        false,                      // flip flux
                        -1,                         // patch for face
                        false,                      // remove from zone
                        -1,                         // zone for face
                        false                       // face zone flip
                    )
                );
            }
            else
            {
                ref.setAction
                (
                    polyModifyFace
                    (
                        newFace,                     // face
                        curFaceID,                   // master face
                        own[curFaceID],              // owner
                        -1,                          // neighbour
                        false,                       // flip flux
                        mesh.boundaryMesh().whichPatch(curFaceID), // patch
                        false,                        // remove from zone
                        -1,                           // zone for face
                        false                         // face zone flip
                    )
                );
            }
        }
    }

    // Break coupled (processor) faces
    forAll (coupledFacesToBreak_, faceI)
    {
        const label& curFaceID = coupledFacesToBreak_[faceI];

        ref.setAction
        (
            polyModifyFace
            (
                faces[curFaceID],            // face
                curFaceID,                   // master face
                own[curFaceID],              // owner
                -1,                          // neighbour
                false,                       // flip flux
                crackPatchID_.index(),       // patch
                false,                       // remove from zone
                -1,                          // zone for face
                false                        // face zone flip
            )
        );
    }

    // Clear the faces to open list
    coupledFacesToBreak_.setSize(0);

    if (debug)
    {
        Pout<< "void faceCracker::detachFaceCracker("
            << "polyTopoChange& ref) const "
            << " for object " << name() << " : "
            << "Finished detaching faces" << endl;
    }
}


// ************************************************************************* //
