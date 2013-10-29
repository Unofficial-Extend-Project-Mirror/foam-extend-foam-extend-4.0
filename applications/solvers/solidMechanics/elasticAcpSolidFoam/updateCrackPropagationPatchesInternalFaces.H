{
    labelHashSet selectedFaces;

    const labelListList& edgeFaces = mesh.edgeFaces();

    const scalarField& cohesiveZoneI = cohesiveZone.internalField();

    forAll (crackPropagationPatches, patchI)
    {
        label curPatchID = crackPropagationPatches[patchI];

        const labelList& patchEdges =
            mesh.boundaryMesh()[curPatchID].meshEdges();

        forAll(patchEdges, edgeI)
        {
            label curEdge = patchEdges[edgeI];

            const labelList& curFaces = edgeFaces[curEdge];

            forAll(curFaces, faceI)
            {
                if (curFaces[faceI] < mesh.nInternalFaces())
                {
                    if (cohesiveZoneI[curFaces[faceI]] > (1.0 - SMALL))
                    {
                        if (!selectedFaces.found(curFaces[faceI]))
                        {
                            selectedFaces.insert(curFaces[faceI]);
                        }
                    }
                }
            }
        }
    }


//     if (notchPatchID != -1)
//     {
//         const labelList& notchPatchEdges =
//             mesh.boundaryMesh()[notchPatchID].meshEdges();

//         forAll(notchPatchEdges, edgeI)
//         {
//             label curEdge = notchPatchEdges[edgeI];

//             const labelList& curFaces = edgeFaces[curEdge];

//             forAll(curFaces, faceI)
//             {
//                 if (curFaces[faceI] < mesh.nInternalFaces())
//                 {
//                     if (cohesiveZoneI[curFaces[faceI]] > (1.0 - SMALL))
//                     {
//                         selectedFaces.insert(curFaces[faceI]);
//                     }
//                 }
//             }
//         }
//     }

//     if (mesh.boundary()[cohesivePatchID].size())
//     {
//         const labelList& cohesivePatchEdges =
//             mesh.boundaryMesh()[cohesivePatchID].meshEdges();

//         forAll(cohesivePatchEdges, edgeI)
//         {
//             label curEdge = cohesivePatchEdges[edgeI];

//             const labelList& curFaces = edgeFaces[curEdge];

//             forAll(curFaces, faceI)
//             {
//                 if (curFaces[faceI] < mesh.nInternalFaces())
//                 {
//                     if (cohesiveZoneI[curFaces[faceI]] > (1.0 - SMALL))
//                     {
//                         selectedFaces.insert(curFaces[faceI]);
//                     }
//                 }
//             }
//         }
//     }

    crackPropagationPatchesInternalFaces = labelList(selectedFaces.toc());
}
