label nCrackTipEdges = 0;

forAll(checkSharedEdge, glEdgeI)
{
    if (checkSharedEdge[glEdgeI])
    {
        label sharedEdgeIndex = findIndex(glEdgeAddr, glEdgeI);

        if ( sharedEdgeIndex != -1 )
//         if ( (sharedEdgeIndex != -1) && (oldCrackTipEdges.size() == 0) )
        {
            label curSharedEdge = glEdges[sharedEdgeIndex];

            const labelList& curFaces = meshEdgeFaces[curSharedEdge];

            label nCrackPatchFaces = 0;

//             vector normal = vector::zero;

            forAll(curFaces, faceI)
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
                oldCrackTipEdges.setSize(nCrackTipEdges+1);
                oldCrackTipEdges[nCrackTipEdges] = curSharedEdge;
            }
        }
    }
}

// Pout << "oldCrackTipEdges_: " << oldCrackTipEdges << endl;
