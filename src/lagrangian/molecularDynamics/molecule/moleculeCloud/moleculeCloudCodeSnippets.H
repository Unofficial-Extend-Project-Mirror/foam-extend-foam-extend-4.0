// Parallel coding to access boundary information to build up interaction cell info

// See preservePatchTypes for how to read the boundary file.

// Read faceProcAddressing, as per reconstructPar, to get hold of the original,
// undecomposed face label from a face on a processor mesh.  See email from Mattijs:

// > Is it a case of reading the faceProcAddressing file, in the same way as
// > something like reconstructPar?
//         Correct.
// 
//         Note that faceProcAddressing is a bit weird since it also includes which side 
// of an internal face we have. If I remember correctly:
// 
//         faceI == 0 illegal
//         faceI > 0 we have the original owner of  faceI-1 i.e. we have the face in the 
//         original order.
//         faceI < 0 we have the original neighbour of -faceI-1 so the face is flipped.

// Use the same functionality as
// label polyBoundaryMesh::whichPatch(const label faceIndex) const
// To determine which patch a face was on originally.

if (Pstream::parRun())
{
//     if (Pstream::myProcNo() == Pstream::masterNo())
// //     {
//         dictionary patchDictionary;
// 
//         DynamicList<word> patchNames;
// 
//         {
//             IOobject undecomposedBoundaryHeader
//             (
//                 "undecomposedBoundary",
//                 mesh_.time().constant(),
//                 polyMesh::meshSubDir,
//                 mesh_,
//                 IOobject::MUST_READ,
//                 IOobject::NO_WRITE,
//                 false
//             );
// 
//             if (undecomposedBoundaryHeader.headerOk())
//             {
//                 polyBoundaryMeshEntries undecomposedPatchEntries
//                 (
//                     undecomposedBoundaryHeader
//                 );
// 
//                 forAll(undecomposedPatchEntries, patchi)
//                 {
//                     patchNames.append
//                     (
//                         undecomposedPatchEntries[patchi].keyword()
//                     );
// 
//                     patchDictionary.add
//                     (
//                         undecomposedPatchEntries[patchi]
//                     );
//                 }
//             }
//             else
//             {
//                 FatalErrorIn
//                 (
//                     "moleculeCloudBuildCellInteractionLists.C\n"
//                 )
//                     << "undecomposedBoundary file not found in "
//                         "constant/polyMesh"
//                     << abort(FatalError);
//             }
//         }
// 
//         labelIOList faceProcAddressing
//         (
//             IOobject
//             (
//                 "faceProcAddressing",
//                 mesh_.time().constant(),
//                 polyMesh::meshSubDir,
//                 mesh_,
//                 IOobject::MUST_READ,
//                 IOobject::NO_WRITE,
//                 false
//             )
//         );

        labelList procPatches(mesh_.globalData().processorPatches());

        forAll(procPatches,pP)
        {
            const processorPolyPatch& patch =
                refCast<const processorPolyPatch>
            (
                mesh_.boundaryMesh()[procPatches[pP]]
            );
// 
//             Pout << nl << "name: " << patch.name() << nl
//                     << "start: " << patch.start() << nl
//                     << "size: " << patch.size() << nl
//                     << "separated: " << Switch(patch.separated()) << nl
//                     << "parallel: " << Switch(patch.parallel()) << nl << endl;
// 
//             forAll (patch, pI)
//             {
//                 label decomposedMeshFace = patch.start() + pI;
// 
//                 label faceProcAdd = faceProcAddressing[decomposedMeshFace];
// 
//                 label globalFace = abs(faceProcAdd)-1;
// 
//                 Pout << "Patch index: " << pI
//                         << " " << patch[pI]
//                         << " Mesh index: " << decomposedMeshFace
//                         << " faceProcAdd: " << faceProcAdd
//                         << " globalFace:" << globalFace;
// 
//                 label minStart = -1;
// 
//                 // Scanning the dictionary each time is a very ugly way of
//                 // finding out what patch a face originally belonged to, but
//                 // it proves the concept.  Read the patch info a container
//                 // class and have a neat way of tell which patch a face is from
//                 // embedded in that.  Split each processor face down into
//                 // separate lists for each different originiating patch.
// 
//                 forAll(patchNames, patchi)
//                 {
//                     if (patchDictionary.found(patchNames[patchi]))
//                     {
//                         const dictionary& patchDict = 
//                                 patchDictionary.subDict(patchNames[patchi]);
// 
//                         word faceName(patchNames[patchi]);
//                         label startFace(readLabel(patchDict.lookup("startFace")));
//                         label nFaces(readLabel(patchDict.lookup("nFaces")));
// 
//                         if
//                         (
//                             minStart < 0
//                          || startFace < minStart
//                         )
//                         {
//                             minStart = startFace;
//                         }
// 
//                         if
//                         (
//                             globalFace >= startFace
//                          && globalFace < startFace + nFaces
//                         )
//                         {
//                             Pout << " original patch: " << faceName  << endl;
//                         }
//                     }
//                 }
// 
//                 if (globalFace < minStart)
//                 {
//                     Pout << " originally an internal face" << endl;
//                 }
//             }
// 
            if (patch.separated())
            {
                 Pout << patch.separation();
            }

            if (!patch.parallel())
            {
                 Pout << patch.forwardT();
            }
        }
//     }
//     else
//     {
// 
//     }

        // Get coords of my shared points
//     vector sharedPoints(vector::one*(Pstream::myProcNo()+1));
//     label testRedLab(Pstream::myProcNo()+1);
    
//     Pout << testRedLab << endl;

        // Append from all processors
//     combineReduce(sharedPoints, plusEqOp<vector>());
//     reduce(testRedLab, plusOp<label>());

//     Pout << testRedLab << endl;
}
