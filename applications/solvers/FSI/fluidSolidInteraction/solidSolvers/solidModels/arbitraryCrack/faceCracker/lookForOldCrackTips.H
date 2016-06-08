label nCrackTipEdges = oldCrackTipEdges.size();
{
    const labelList& curFaceEdges = mesh.faceEdges(curFaceID);

//     const label curFaceCell = own[curFaceID];

//     Pout << "curFaceEdges: " << curFaceEdges << endl;

    forAll(curFaceEdges, edgeI)
    {
        label curEdge = curFaceEdges[edgeI];

        const labelList& curFaces = meshEdgeFaces[curEdge];

        label nCrackPatchFaces = 0;

//         vector normal = vector::zero;

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
            oldCrackTipEdges[nCrackTipEdges] = curEdge;
            nCrackTipEdges++;
        }
        else
        {
            label glEdge =findIndex(glEdges, curEdge);

            if (glEdge != -1)
            {
                checkSharedEdge[glEdgeAddr[glEdge]] = 1;
            }
        }
    }
}

// Pout << "oldCrackTipEdges: " << oldCrackTipEdges << endl;
