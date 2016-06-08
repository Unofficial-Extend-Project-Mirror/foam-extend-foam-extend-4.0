const edgeList& edges = mesh.edges();

forAll(oldCrackTipEdges, edgeI)
{
    label curEdge = oldCrackTipEdges[edgeI];

    label masterCrackFace = -1;
    label slaveCrackFace = -1;

    const labelList& curFaces = meshEdgeFaces[curEdge];

    forAll(curFaces, faceI)
    {
        if (!mesh.isInternalFace(curFaces[faceI]))
        {
            label patchID =
                mesh.boundaryMesh().whichPatch(curFaces[faceI]);

//             label localFaceID =
//                 curFaces[faceI] - mesh.boundaryMesh()[patchID].start()

            if (patchID == crackPatchID_.index())
            {
//                 if
//                 (
//                     (
//                         oldCrackTipMasterFaceNormal[edgeI]
//                       & n[curFaces[faceI]]
//                     ) > SMALL
//                 )
//                 {
//                     masterCrackFace = curFaces[faceI];
//                 }
//                 else
//                 {
//                     slaveCrackFace = curFaces[faceI];
//                 }

                if (masterCrackFace == -1)
                {
                    masterCrackFace = curFaces[faceI];
                }
                else
                {
                    slaveCrackFace = curFaces[faceI];
                }
            }
        }
    }

    label masterCell = own[masterCrackFace];
    label slaveCell = own[slaveCrackFace];

    // Create new points
    labelList addedPoints(2, -1);

    addedPoints[0] =
        ref.setAction
        (
            polyAddPoint
            (
                points[edges[curEdge].start()],        // point
                edges[curEdge].start(),                // master point
                -1,                        // zone ID
                true                       // supports a cell
            )
        );

    addedPoints[1] =
        ref.setAction
        (
            polyAddPoint
            (
                points[edges[curEdge].end()],        // point
                edges[curEdge].end(),                // master point
                -1,                        // zone ID
                true                       // supports a cell
            )
        );

//     Pout << "added crack tip points: " << addedPoints << endl;

    labelHashSet masterCellFaceMap(6);

    // First we add all pointFaces of new points
    const labelListList& pointFaces = mesh.pointFaces();

    forAll(pointFaces[edges[curEdge].start()], fI)
    {
        label pointFaceI = pointFaces[edges[curEdge].start()][fI];
        masterCellFaceMap.insert(pointFaceI);
    }

    forAll(pointFaces[edges[curEdge].end()], fI)
    {
        label pointFaceI = pointFaces[edges[curEdge].end()][fI];
        masterCellFaceMap.insert(pointFaceI);
    }

    SLList<label> cellsToCheck;
    cellsToCheck.append(slaveCell);

    do
    {
        // We will remove all faces of the current cell from the
        // masterCellFaceMap if they are found, and add the neighbour cell
        // the cellsToCheck list if it is not the master cell
        const label cellID = cellsToCheck.first();
        const labelList& curFaces = cells[cellID];
        forAll(curFaces, faceI)
        {
            const label curFaceI = curFaces[faceI];
            if (masterCellFaceMap.found(curFaceI))
            {
                masterCellFaceMap.erase(curFaceI);

                // Add neighbour cell to cellsToCheck if it is not the
                // master cell
                if (mesh.isInternalFace(curFaceI))
                {
                    label neiCell = nei[curFaceI];
                    if (neiCell == cellID)
                    {
                        neiCell = own[curFaceI];
                    }
                    bool neiIsMaster = false;
                    if (masterCell == neiCell)
                    {
                        neiIsMaster = true;
                    }

                    if (!neiIsMaster)
                    {
                        cellsToCheck.append(neiCell);
                    }
                }
            }
        }

        // Remove current cell from cellsToCheck
        cellsToCheck.removeHead();
    }
    while (cellsToCheck.size());


    // Create the master layer point map
    Map<label> masterLayerPointMap(2);

    masterLayerPointMap.insert
    (
        edges[curEdge].start(),
        addedPoints[0]
    );
    masterLayerPointMap.insert
    (
        edges[curEdge].end(),
        addedPoints[1]
    );


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
            // We loop through all points of the masterCellFace and if a point
            // from the cracked face is found then this face needs to be
            // changed

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
                label coupledFaceToBreakIndex =
                    findIndex(coupledFacesToBreak_, curFaceID);

                if (coupledFaceToBreakIndex == -1)
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
                            crackPatchID_.index(),       // patch
                            false,                       // remove from zone
                            -1,                          // zone for face
                            false                        // face zone flip
                        )
                    );
                }
            }
        }
    }

    forAll (coupledFacesToBreak_, faceI)
    {
        const label& curFaceID = coupledFacesToBreak_[faceI];

        if (!masterCellFaceMap.found(curFaceID))
        {
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
    }
}
