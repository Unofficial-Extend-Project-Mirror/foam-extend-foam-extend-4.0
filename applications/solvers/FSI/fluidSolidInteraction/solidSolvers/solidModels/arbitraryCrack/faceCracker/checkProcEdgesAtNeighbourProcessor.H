{
    List<labelPair> procPatchEdge(procPatchEdgeSet);
    List<labelPair> receivedProcPatchEdge(0);

    // Old crack tip edges at the processor boundary
    // whose crack patch faces should be detached
    // because new crack tip is now on the neighbour processor
    label nCrackTipEdges = oldCrackTipEdges.size();

    // Send edges
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<processorPolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh.boundaryMesh()[patchI]
                );

            labelHashSet curPatchEdgesSet;

            forAll (procPatchEdge, edgeI)
            {
                if (procPatchEdge[edgeI].first() == patchI)
                {
                    if
                    (
                       !curPatchEdgesSet.found
                        (
                            procPatchEdge[edgeI].second()
                        )
                    )
                    {
                        curPatchEdgesSet.insert
                        (
                            procPatchEdge[edgeI].second()
                        );
                    }
                }
            }

            labelList curPatchEdges = curPatchEdgesSet.toc();

            OPstream toProc
            (
                Pstream::blocking,
                procPatch.neighbProcNo()
            );

            toProc << curPatchEdges;
        }
    }

    // Receive edges
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<processorPolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh.boundaryMesh()[patchI]
                );

            labelList ngbProcEdges;

            IPstream fromProc
            (
                Pstream::blocking,
                procPatch.neighbProcNo()
            );

            fromProc >> ngbProcEdges;

            forAll(ngbProcEdges, edgeI)
            {
                label curSize = receivedProcPatchEdge.size();
                receivedProcPatchEdge.setSize(curSize+1);
                receivedProcPatchEdge[curSize] =
                    labelPair
                    (
                        patchI,
                        findIndex
                        (
                            procPatch.neighbEdges(),
                            ngbProcEdges[edgeI]
                        )
                    );
            }
        }
    }

    List<bool> receivedProcPatchEdgeIsInternal
    (
        receivedProcPatchEdge.size(),
        true
    );

    forAll(receivedProcPatchEdge, edgeI)
    {
        label curPatchID = receivedProcPatchEdge[edgeI].first();
        label curEdgeID = receivedProcPatchEdge[edgeI].second();

//         receivedProcPatchEdgeIsInternal[edgeI] = true;

        const labelList& curMeshEdges =
            mesh.boundaryMesh()[curPatchID].meshEdges();

        const labelList& curFaces =
            meshEdgeFaces[curMeshEdges[curEdgeID]];

        forAll (curFaces, faceI)
        {
            if (!mesh.isInternalFace(curFaces[faceI]))
            {
                label patchID =
                    mesh.boundaryMesh().whichPatch(curFaces[faceI]);

                if (!isA<processorPolyPatch>(mesh.boundaryMesh()[patchID]))
                {
                    // The edge belongs to a boundary face
                    receivedProcPatchEdgeIsInternal[edgeI] = false;
                    break;
                }
            }
        }

        if (!receivedProcPatchEdgeIsInternal[edgeI])
        {
            label nCrackPatchFaces = 0;

            forAll (curFaces, faceI)
            {
                if (!mesh.isInternalFace(curFaces[faceI]))
                {
                    label patchID =
                        mesh.boundaryMesh().whichPatch(curFaces[faceI]);

                    if (patchID == crackPatchID_.index())
                    {
                        nCrackPatchFaces++;
                    }
                }
            }

            if (nCrackPatchFaces == 2)
            {
                oldCrackTipEdges.setSize(nCrackTipEdges + 1);

                oldCrackTipEdges[nCrackTipEdges] = curMeshEdges[curEdgeID];

                nCrackTipEdges++;
            }
        }
    }

//     oldCrackTipEdges.setSize(nCrackTipEdges);

    // Send edges
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<processorPolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh.boundaryMesh()[patchI]
                );

            labelList curPatchEdges
            (
                receivedProcPatchEdgeIsInternal.size(),
                -1
            );
            boolList curPatchEdgeIsInternal
            (
                receivedProcPatchEdgeIsInternal.size(),
                false
            );

            label index = 0;
            forAll(receivedProcPatchEdgeIsInternal, edgeI)
            {
                if (receivedProcPatchEdge[edgeI].first() == patchI)
                {
                    curPatchEdges[index] =
                        receivedProcPatchEdge[edgeI].second();
                    curPatchEdgeIsInternal[index] =
                        receivedProcPatchEdgeIsInternal[edgeI];
                    index++;
                }
            }
            curPatchEdges.setSize(index);
            curPatchEdgeIsInternal.setSize(index);

            OPstream toProc
            (
                Pstream::blocking,
                procPatch.neighbProcNo()
            );

            toProc << curPatchEdges << curPatchEdgeIsInternal << endl;
        }
    }


    // Receive edges and check if internal
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<processorPolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh.boundaryMesh()[patchI]
                );

            labelList ngbProcEdges;
            boolList ngbProcEdgeIsInternal;

            IPstream fromProc
            (
                Pstream::blocking,
                procPatch.neighbProcNo()
            );

            fromProc >> ngbProcEdges >> ngbProcEdgeIsInternal;

            forAll(ngbProcEdges, edgeI)
            {
                label curEdge =
                    findIndex
                    (
                        procPatch.neighbEdges(),
                        ngbProcEdges[edgeI]
                    );

                bool isInternal = ngbProcEdgeIsInternal[edgeI];

                if (isInternal)
                {
                    label curEdgeID =
                        findIndex
                        (
                            meshEdges,
                            procPatch.meshEdges()[curEdge]
                        );

                    // Reset the point creation
                    addedPoints[zoneLocalEdges[curEdgeID].start()] =
                        mp[zoneLocalEdges[curEdgeID].start()];

                    addedPoints[zoneLocalEdges[curEdgeID].end()] =
                        mp[zoneLocalEdges[curEdgeID].end()];
                }
            }
        }
    }
}

{
    List<label> sharedEdges(sharedEdgeSet);

    label nCrackTipEdges = oldCrackTipEdges.size();

    labelList checkSharedEdge(mesh.globalData().nGlobalEdges(), 0);

    forAll(sharedEdges, edgeI)
    {
        label curEdge = sharedEdges[edgeI];

        label glEdge = findIndex(glEdges, curEdge);

        if (glEdge != -1)
        {
            checkSharedEdge[glEdgeAddr[glEdge]] = 1;
        }
    }

    reduce(checkSharedEdge, maxOp<labelList>());

    labelList internalSharedEdge(mesh.globalData().nGlobalEdges(), 1);

    forAll(checkSharedEdge, glEdgeI)
    {
        if (checkSharedEdge[glEdgeI])
        {
            label sharedEdgeIndex = findIndex(glEdgeAddr, glEdgeI);

            if ( sharedEdgeIndex != -1 )
            {
                label curSharedEdge = glEdges[sharedEdgeIndex];

                const labelList& curFaces = meshEdgeFaces[curSharedEdge];

                label nCrackPatchFaces = 0;

                forAll(curFaces, faceI)
                {
                    if (!mesh.isInternalFace(curFaces[faceI]))
                    {
                        label patchID =
                            mesh.boundaryMesh().whichPatch(curFaces[faceI]);

                        if
                        (
                           !isA<processorPolyPatch>
                            (
                                mesh.boundaryMesh()[patchID]
                            )
                        )
                        {
                            // The edge belongs to a boundary face
                            internalSharedEdge[glEdgeI] = 0;
                        }

                        if (patchID == crackPatchID_.index())
                        {
                            nCrackPatchFaces++;
                        }
                    }
                }

                if (nCrackPatchFaces == 2)
                {
                    oldCrackTipEdges.setSize(nCrackTipEdges+1);

                    oldCrackTipEdges[nCrackTipEdges] = curSharedEdge;
                }
            }
        }
    }

    reduce(internalSharedEdge, minOp<labelList>());

    forAll(sharedEdges, edgeI)
    {
        label curEdge = sharedEdges[edgeI];

        label glEdge = findIndex(glEdges, curEdge);

        label isInternal = 1;

        if (glEdge != -1)
        {
            isInternal = internalSharedEdge[glEdgeAddr[glEdge]];
        }

//         Pout << "isInternal " << isInternal << endl;

        if (isInternal)
        {
            label curEdgeID =
                findIndex
                (
                    meshEdges,
                    curEdge
                );

            // Reset the point creation
            addedPoints[zoneLocalEdges[curEdgeID].start()] =
                mp[zoneLocalEdges[curEdgeID].start()];

            addedPoints[zoneLocalEdges[curEdgeID].end()] =
                mp[zoneLocalEdges[curEdgeID].end()];
        }
    }
}

// Pout << "Old crack tip edges xxxxx: " << oldCrackTipEdges << endl;
