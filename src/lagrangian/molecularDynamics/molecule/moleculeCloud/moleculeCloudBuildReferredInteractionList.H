Info << nl << "Building list of referred interaction neighbours" << endl;

labelListList processorPatchSegmentMapping
(
    mesh_.globalData().processorPatches().size()
);

List<vectorList> allNeighbourFaceCentres
(
    mesh_.globalData().processorPatches().size()
);

List<vectorList> allNeighbourFaceAreas
(
    mesh_.globalData().processorPatches().size()
);

label nUndecomposedPatches = 0;

if (Pstream::parRun())
{
    dictionary patchDictionary;

    DynamicList<word> patchNames;

//     IOobject undecomposedBoundaryHeader
//     (
//         "undecomposedBoundary",
//         mesh_.time().constant(),
//         polyMesh::meshSubDir,
//         mesh_,
//         IOobject::MUST_READ,
//         IOobject::NO_WRITE,
//         false
//     );

    Time undecomposedTime
    (
        Time::controlDictName,
        mesh_.time().rootPath(),
        mesh_.time().caseName().path()
    );

    IOobject undecomposedBoundaryHeader
    (
        "boundary",
        undecomposedTime.constant(),
        polyMesh::meshSubDir,
        undecomposedTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (undecomposedBoundaryHeader.headerOk())
    {
        polyBoundaryMeshEntries undecomposedPatchEntries
        (
            undecomposedBoundaryHeader
        );

        forAll(undecomposedPatchEntries, patchi)
        {
            patchNames.append
            (
                undecomposedPatchEntries[patchi].keyword()
            );

            patchDictionary.add
            (
                undecomposedPatchEntries[patchi]
            );
        }
    }
    else
    {
        FatalErrorIn
        (
            "moleculeCloudBuildCellInteractionLists.C"
        )
            << nl << "unable to read undecomposed boundary file from "
            "constant/polyMesh" << nl
            << abort(FatalError);
    }

    labelIOList faceProcAddressing
    (
        IOobject
        (
            "faceProcAddressing",
            mesh_.time().constant(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    labelList procPatches(mesh_.globalData().processorPatches());

    nUndecomposedPatches = patchNames.size();

    // processorPatchSegmentMapping works by mapping the original patch
    // and half that a face on a processor patch was on before decomposition.
    // This creates a patch segment for each half of each original (cyclic)
    // patch which can be assessed separately.
    // There are n = patchNames.size() original patches, k = 0 to n-1.
    // The mapping is:
    // value = 0:  originally an internal face
    // value = k, was originally on the on the patch k-1, 1st half
    // value = -k, was originally on the on the patch k-1, 2nd half

    forAll(procPatches,pP)
    {
        const processorPolyPatch& patch =
            refCast<const processorPolyPatch>
            (
                mesh_.boundaryMesh()[procPatches[pP]]
            );

        labelList& procPatchSegMap = processorPatchSegmentMapping[pP];

        procPatchSegMap.setSize(patch.size());

        forAll (patch, pI)
        {
            label decomposedMeshFace = patch.start() + pI;

            label faceProcAdd = faceProcAddressing[decomposedMeshFace];

            label globalFace = abs(faceProcAdd)-1;

//             Pout << "Patch index: " << pI
//                     << " " << patch[pI]
//                     << " Mesh index: " << decomposedMeshFace
//                     << " faceProcAdd: " << faceProcAdd
//                     << " globalFace:" << globalFace << endl;

            label minStart = -1;

            forAll(patchNames, pN)
            {
                if (patchDictionary.found(patchNames[pN]))
                {
                    const dictionary& patchDict =
                        patchDictionary.subDict(patchNames[pN]);

//                     word faceName(patchNames[pN]);
                    label startFace(readLabel(patchDict.lookup("startFace")));
                    label nFaces(readLabel(patchDict.lookup("nFaces")));

                    if (minStart < 0 || startFace < minStart)
                    {
                        minStart = startFace;
                    }

                    if
                    (
                        globalFace >= startFace
                     && globalFace < startFace + nFaces/2
                    )
                    {
                        procPatchSegMap[pI] = pN + 1;
                    }
                    else if
                    (
                        globalFace >= startFace + nFaces/2
                     && globalFace < startFace + nFaces
                    )
                    {
                        procPatchSegMap[pI] = -(pN + 1);
                    }
                }
            }

            if (globalFace < minStart)
            {
                procPatchSegMap[pI] = 0;
            }
        }
    }

    forAll(procPatches,pP)
    {
        const processorPolyPatch& patch =
            refCast<const processorPolyPatch>
            (
                mesh_.boundaryMesh()[procPatches[pP]]
            );

        {
            OPstream toNeighbProc
            (
                Pstream::blocking,
                patch.neighbProcNo(),
                + 2*(sizeof(label) + patch.size()*sizeof(vector))
            );

            toNeighbProc
                << patch.faceCentres()
                << patch.faceAreas();

//             Pout << "sending face info to processor "
//                 << patch.neighbProcNo() << ", size = "
//                 << patch.size() << endl;
        }
    }

    forAll(procPatches,pP)
    {
        const processorPolyPatch& patch =
            refCast<const processorPolyPatch>
            (
                mesh_.boundaryMesh()[procPatches[pP]]
            );

        vectorList& neighbFaceCentres = allNeighbourFaceCentres[pP];

        neighbFaceCentres.setSize(patch.size());

        vectorList& neighbFaceAreas = allNeighbourFaceAreas[pP];

        neighbFaceAreas.setSize(patch.size());

        {
            IPstream fromNeighbProc
            (
                Pstream::blocking,
                patch.neighbProcNo(),
                2*(sizeof(label) + patch.size()*sizeof(vector) + sizeof(scalar))
            );
            fromNeighbProc
                >> neighbFaceCentres
                >> neighbFaceAreas;

//             Pout << "receiving face info from processor "
//                 << patch.neighbProcNo() << ", size = "
//                 << patch.size() << endl;
        }
    }

    // *************************************************************
    // Tests that all processor patch segments from different
    // original patches prior to decomposition produce the same
    // transform. Check before 1st iteration.
    // *************************************************************

    forAll(procPatches,pP)
    {
        const processorPolyPatch& patch =
            refCast<const processorPolyPatch>
            (
                mesh_.boundaryMesh()[procPatches[pP]]
            );

        const vectorList& neighbFaceCentres = allNeighbourFaceCentres[pP];

        const vectorList& neighbFaceAreas = allNeighbourFaceAreas[pP];

        label nUP;

        for
        (
            nUP = -nUndecomposedPatches;
            nUP <= nUndecomposedPatches;
            nUP++
        )
        {
            DynamicList<vector> refOff;

            DynamicList<tensor> refTrans;

            forAll (patch, faceI)
            {
                if (processorPatchSegmentMapping[pP][faceI] == nUP)
                {
                    referredCell testRefCell
                    (
                        mesh_,
                        -1,
                        -1,
                        patch.faceCentres()[faceI],
                        neighbFaceCentres[faceI],
                        patch.faceNormals()[faceI],
                        neighbFaceAreas[faceI]
                       /(mag(neighbFaceAreas[faceI]) + VSMALL)
                    );

                    refOff.append(testRefCell.offset());

                    refTrans.append(testRefCell.rotation());
                }
            }

            refOff.shrink();

            refTrans.shrink();

            if (refOff.size())
            {
                if
                (
                    sum(mag(refOff-refOff[0]))/refOff.size() > transTol
                 || sum(mag(refTrans-refTrans[0]))/refTrans.size() > transTol
                )
                {
                    FatalErrorIn
                    (
                        "moleculeCloud::buildCellInteractionLists"
                    )
                        << nl << "Face pairs on patch "
                        << patch.name()
                        << ", segment " << patchNames[nUP]
                        << " do not give the same referring "
                        << " transformations to within tolerance of "
                        << transTol << nl
                        << " Referring offsets:" << refOff << nl
                        << " Average sum of mag difference: "
                        << sum(mag(refOff-refOff[0]))/refOff.size() << nl
                        << " Referring transforms:" << refTrans << nl
                        << " Average sum of mag difference: "
                        << sum(mag(refTrans-refTrans[0]))/refTrans.size()
                        << nl << abort(FatalError);
                }
            }
        }
    }
}

label cellsReferredThisIteration = 1;

label iterationNo = 0;

DynamicList<label> realCellsWithinRCutMaxOfAnyReferringPatch;

while(cellsReferredThisIteration)
{
    label refIntListStartSize = referredInteractionList.size();

    forAll (mesh_.boundaryMesh(), patchI)
    {
        // Treat local cyclics on each processor before processor boundaries.
        // Separate treatment allows the serial version to run transparently.

        if (mesh_.boundaryMesh()[patchI].type() == "cyclic")
        {
            const cyclicPolyPatch& patch =
                    refCast<const cyclicPolyPatch>
            (
                mesh_.boundaryMesh()[patchI]
            );

            if(patch.size())
            {
    //             Pout << patch << endl;

                if (iterationNo == 0)
                {
                    // Tests that all combinations of face pairs produce the
                    // same transform.  Only check on 1st iteration.

                    label faceL;
                    // A face in the 1st half of the patch

                    label faceM;
                    // Face corresponding to faceL in the 2nd half of the
                    // patch. Assumes correct face ordering.

                    vectorList refOff(patch.size()/2);

                    List<tensor> refTrans(patch.size()/2);

                    for
                    (
                        faceL = 0, faceM = patch.size()/2;
                        faceL < patch.size()/2;
                        faceL++, faceM++
                    )
                    {
                        referredCell testRefCell
                        (
                            mesh_,
                            -1,
                            -1,
                            patch.faceCentres()[faceL],
                            patch.faceCentres()[faceM],
                            patch.faceNormals()[faceL],
                            patch.faceNormals()[faceM]
                        );

                        refOff[faceL] = testRefCell.offset();

                        refTrans[faceL] = testRefCell.rotation();
                    }

                    if
                    (
                        sum(mag(refOff - refOff[0]))/(patch.size()/2)
                      > transTol
                     || sum(mag(refTrans - refTrans[0]))/(patch.size()/2)
                      > transTol
                    )
                    {
                        FatalErrorIn
                        (
                            "moleculeCloud::buildCellInteractionLists"
                        )
                            << nl << "Face pairs on patch "
                            << patch.name()
                            << " do not give the same referring "
                            << " transformations to within tolerance of "
                            << transTol << nl
                            << " Referring offsets:" << refOff << nl
                            << " Average sum of mag difference: "
                            << sum(mag(refOff-refOff[0]))/refOff.size() << nl
                            << " Referring transforms:" << refTrans << nl
                            << " Average sum of mag difference: "
                            << sum(mag(refTrans-refTrans[0]))/refTrans.size()
                            << nl << abort(FatalError);
                    }
                }

                // ***********************************************************
                // 1st half of face list - 1st side of boundary
                // ***********************************************************

                label faceI;

                DynamicList<label> meshFacesOnThisSegment;

                for (faceI = 0; faceI < patch.size()/2; faceI++)
                {
                    // unable to use the normal accessors for the polyPatch
                    // because points on separate halves need used separately

                    meshFacesOnThisSegment.append(faceI + patch.start());
                }

                meshFacesOnThisSegment.shrink();

                DynamicList<label> meshEdgesOnThisSegment;

                DynamicList<label> meshPointsOnThisSegment;

                forAll(meshFacesOnThisSegment, mFOTS)
                {
                    const label segFace = meshFacesOnThisSegment[mFOTS];

                    const labelList& faceEdges = mesh_.faceEdges()[segFace];

                    forAll (faceEdges, fE)
                    {
                        const label faceEdge(faceEdges[fE]);

                        if (findIndex(meshEdgesOnThisSegment, faceEdge) == -1)
                        {
                            meshEdgesOnThisSegment.append(faceEdge);
                        }
                    }

                    const face& facePoints(mesh_.faces()[segFace]);

                    forAll (facePoints, fP)
                    {
                        const label facePoint(facePoints[fP]);

                        if
                        (
                            findIndex(meshPointsOnThisSegment, facePoint)
                         == -1
                        )
                        {
                            meshPointsOnThisSegment.append(facePoint);
                        }
                    }
                }

                meshEdgesOnThisSegment.shrink();

                meshPointsOnThisSegment.shrink();

                if (iterationNo == 0)
                {
                    // Assessing real cells in range is only required on
                    // the 1st iteration because they do not change from
                    // iteration to iteration.

                    labelList realCellsFoundInRange
                    (
                        realCellsInRangeOfSegment
                        (
                            meshFacesOnThisSegment,
                            meshEdgesOnThisSegment,
                            meshPointsOnThisSegment
                        )
                    );

                    forAll(realCellsFoundInRange,cFIR)
                    {
                        const label realCell = realCellsFoundInRange[cFIR];

                        referredCell cellToRefer
                        (
                            mesh_,
                            Pstream::myProcNo(),
                            realCell,
                            patch.faceCentres()[0],
                            patch.faceCentres()[patch.size()/2],
                            patch.faceNormals()[0],
                            patch.faceNormals()[patch.size()/2]
                        );

                        // Test all existing referred and real cells to
                        // check duplicates are not being made or cells
                        // aren't being referred back onto themselves

                        bool addCellToRefer = true;

                        // Check if cellToRefer is an existing referred cell

                        forAll(referredInteractionList, rIL)
                        {
                            if
                            (
                                cellToRefer.duplicate
                                (
                                    referredInteractionList[rIL]
                                )
                            )
                            {
                                addCellToRefer = false;

                                break;
                            }
                        }

                        // Check for cellToRefer being referred back
                        // ontop of a real cell

                        if
                        (
                            cellToRefer.duplicate
                            (
                                Pstream::myProcNo(),
                                mesh_.nCells()
                            )
                        )
                        {
                            addCellToRefer = false;
                        }

                        if (addCellToRefer)
                        {
                            referredInteractionList.append(cellToRefer);
                        }

                        // add real cells found in range of cyclic patch
                        // to whole mesh list

                        if
                        (
                            findIndex
                            (
                                realCellsWithinRCutMaxOfAnyReferringPatch,
                                realCell
                            ) == -1
                        )
                        {
                            realCellsWithinRCutMaxOfAnyReferringPatch.append
                            (
                                realCell
                            );
                        }
                    }
                }

                referredInteractionList.shrink();

                labelList referredCellsFoundInRange
                (
                    referredCellsInRangeOfSegment
                    (
                        referredInteractionList,
                        meshFacesOnThisSegment,
                        meshEdgesOnThisSegment,
                        meshPointsOnThisSegment
                    )
                );

                forAll(referredCellsFoundInRange,cFIR)
                {
                    referredCell& existingRefCell = referredInteractionList
                    [
                        referredCellsFoundInRange[cFIR]
                    ];

                    referredCell cellToReRefer = existingRefCell.reRefer
                    (
                        patch.faceCentres()[0],
                        patch.faceCentres()[patch.size()/2],
                        patch.faceNormals()[0],
                        patch.faceNormals()[patch.size()/2]
                    );

                    // Test all existing referred and real cells to check
                    // duplicates are not being made or cells aren't being
                    // referred back onto themselves

                    bool addCellToReRefer = true;

                    // Check if cellToRefer is an existing referred cell

                    forAll(referredInteractionList, rIL)
                    {
                        if
                        (
                            cellToReRefer.duplicate
                            (
                                referredInteractionList[rIL]
                            )
                        )
                        {
                            addCellToReRefer = false;

                            break;
                        }
                    }

                    // Check for cellToRefer being referred back
                    // ontop of a real cell

                    if
                    (
                        cellToReRefer.duplicate
                        (
                            Pstream::myProcNo(),
                            mesh_.nCells()
                        )
                    )
                    {
                        addCellToReRefer = false;
                    }

                    if(addCellToReRefer)
                    {
                        referredInteractionList.append(cellToReRefer);
                    }
                }

                // ***********************************************************
                // 2nd half of face list - 2nd side of boundary
                // ***********************************************************

                meshFacesOnThisSegment.clear();

                for (faceI = patch.size()/2; faceI < patch.size(); faceI++)
                {
                    // unable to use the normal accessors for the polyPatch
                    // because points on separate halves need used separately

                    meshFacesOnThisSegment.append(faceI + patch.start());
                }

                meshFacesOnThisSegment.shrink();

                meshEdgesOnThisSegment.clear();

                meshPointsOnThisSegment.clear();

                forAll(meshFacesOnThisSegment, mFOTS)
                {
                    const label segFace = meshFacesOnThisSegment[mFOTS];

                    const labelList& faceEdges = mesh_.faceEdges()[segFace];

                    forAll (faceEdges, fE)
                    {
                        const label faceEdge(faceEdges[fE]);

                        if(findIndex(meshEdgesOnThisSegment, faceEdge) == -1)
                        {
                            meshEdgesOnThisSegment.append(faceEdge);
                        }
                    }

                    const face& facePoints(mesh_.faces()[segFace]);

                    forAll (facePoints, fP)
                    {
                        const label facePoint(facePoints[fP]);

                        if
                        (
                           findIndex(meshPointsOnThisSegment, facePoint)
                         == -1
                        )
                        {
                            meshPointsOnThisSegment.append(facePoint);
                        }
                    }
                }

                meshEdgesOnThisSegment.shrink();

                meshPointsOnThisSegment.shrink();

                if (iterationNo == 0)
                {
                    // Assessing real cells in range is only required on
                    // the 1st iteration because they do not change from
                    // iteration to iteration.

                    labelList realCellsFoundInRange
                    (
                        realCellsInRangeOfSegment
                        (
                            meshFacesOnThisSegment,
                            meshEdgesOnThisSegment,
                            meshPointsOnThisSegment
                        )
                    );

                    forAll(realCellsFoundInRange,cFIR)
                    {
                        const label realCell = realCellsFoundInRange[cFIR];

                        referredCell cellToRefer
                        (
                            mesh_,
                            Pstream::myProcNo(),
                            realCell,
                            patch.faceCentres()[patch.size()/2],
                            patch.faceCentres()[0],
                            patch.faceNormals()[patch.size()/2],
                            patch.faceNormals()[0]
                        );

                        // Test all existing referred and real cells to
                        // check duplicates are not being made or cells
                        // aren't being referred back onto themselves

                        bool addCellToRefer = true;

                        // Check if cellToRefer is an existing referred cell

                        forAll(referredInteractionList, rIL)
                        {
                            if
                            (
                                cellToRefer.duplicate
                                (
                                    referredInteractionList[rIL]
                                )
                            )
                            {
                                addCellToRefer = false;

                                break;
                            }
                        }

                        // Check for cellToRefer being referred back
                        // ontop of a real cell

                        if
                        (
                            cellToRefer.duplicate
                            (
                                Pstream::myProcNo(),
                                mesh_.nCells()
                            )
                        )
                        {
                            addCellToRefer = false;
                        }

                        if (addCellToRefer)
                        {
                            referredInteractionList.append(cellToRefer);
                        }

                        // add real cells found in range of cyclic patch
                        // to whole mesh list

                        if
                        (
                            findIndex
                            (
                                realCellsWithinRCutMaxOfAnyReferringPatch,
                                realCell
                            ) == -1
                        )
                        {
                            realCellsWithinRCutMaxOfAnyReferringPatch.append
                            (
                                realCell
                            );
                        }
                    }
                }

                referredInteractionList.shrink();

                referredCellsFoundInRange = referredCellsInRangeOfSegment
                (
                    referredInteractionList,
                    meshFacesOnThisSegment,
                    meshEdgesOnThisSegment,
                    meshPointsOnThisSegment
                );

                forAll(referredCellsFoundInRange,cFIR)
                {
                    referredCell& existingRefCell = referredInteractionList
                    [
                        referredCellsFoundInRange[cFIR]
                    ];

                    referredCell cellToReRefer = existingRefCell.reRefer
                    (
                        patch.faceCentres()[patch.size()/2],
                        patch.faceCentres()[0],
                        patch.faceNormals()[patch.size()/2],
                        patch.faceNormals()[0]
                    );

                    // Test all existing referred and real cells to check
                    // duplicates are not being made or cells aren't being
                    // referred back onto themselves

                    bool addCellToReRefer = true;

                    // Check if cellToRefer is an existing referred cell

                    forAll(referredInteractionList, rIL)
                    {
                        if
                        (
                            cellToReRefer.duplicate
                            (
                                referredInteractionList[rIL]
                            )
                        )
                        {
                            addCellToReRefer = false;

                            break;
                        }
                    }

                    // Check for cellToRefer being referred back
                    // ontop of a real cell

                    if
                    (
                        cellToReRefer.duplicate
                        (
                            Pstream::myProcNo(),
                            mesh_.nCells()
                        )
                    )
                    {
                        addCellToReRefer = false;
                    }

                    if (addCellToReRefer)
                    {
                        referredInteractionList.append(cellToReRefer);
                    }
                }
            }
        }
    }

    if (Pstream::parRun())
    {
        labelList procPatches(mesh_.globalData().processorPatches());

        forAll(procPatches,pP)
        {
            const processorPolyPatch& patch =
                refCast<const processorPolyPatch>
                (
                    mesh_.boundaryMesh()[procPatches[pP]]
                );

            DynamicList<referredCell> referredCellsToTransfer;

            const vectorList& neighbFaceCentres = allNeighbourFaceCentres[pP];

            const vectorList& neighbFaceAreas = allNeighbourFaceAreas[pP];

            label nUP;

            for
            (
                nUP = -nUndecomposedPatches;
                nUP <= nUndecomposedPatches;
                nUP++
            )
            {
                // faceT is used to specify one face on this patch segment that
                // will be used to calculate the transformation values.  All
                // faces are guaranteed to produce the same transform because
                // of the checks carried out at the start of the function.
                // Setting to -1 until the 1st face on this segment is found.

                label faceT = -1;

                DynamicList<label> meshFacesOnThisSegment;

                forAll (patch, faceI)
                {
                    if (processorPatchSegmentMapping[pP][faceI] == nUP)
                    {
                        if (faceT == -1)
                        {
                            faceT = faceI;
                        }

                        meshFacesOnThisSegment.append(faceI + patch.start());
                    }
                }

                meshFacesOnThisSegment.shrink();

                DynamicList<label> meshEdgesOnThisSegment;

                DynamicList<label> meshPointsOnThisSegment;

                forAll(meshFacesOnThisSegment, mFOTS)
                {
                    const label segFace = meshFacesOnThisSegment[mFOTS];

                    const labelList& faceEdges = mesh_.faceEdges()[segFace];

                    forAll (faceEdges, fE)
                    {
                        const label faceEdge(faceEdges[fE]);

                        if(findIndex(meshEdgesOnThisSegment, faceEdge) == -1)
                        {
                            meshEdgesOnThisSegment.append(faceEdge);
                        }
                    }

                    const face& facePoints(mesh_.faces()[segFace]);

                    forAll (facePoints, fP)
                    {
                        const label facePoint(facePoints[fP]);

                        if
                        (
                           findIndex(meshPointsOnThisSegment, facePoint)
                         == -1
                        )
                        {
                            meshPointsOnThisSegment.append(facePoint);
                        }
                    }
                }

                meshEdgesOnThisSegment.shrink();

                meshPointsOnThisSegment.shrink();

                if (meshFacesOnThisSegment.size())
                {
                    if (faceT == -1)
                    {
                        FatalErrorIn
                        (
                            "moleculeCloudBuildReferredInteractionList.H"
                        )
                        << nl << "faceT == -1 encountered but "
                        << meshFacesOnThisSegment.size()
                        << " faces found on patch segment."
                        << abort(FatalError);
                    }

                    if (iterationNo == 0)
                    {
                        // Assessing real cells in range is only required on
                        // the 1st iteration because they do not change from
                        // iteration to iteration.

                        labelList realCellsFoundInRange
                        (
                            realCellsInRangeOfSegment
                            (
                                meshFacesOnThisSegment,
                                meshEdgesOnThisSegment,
                                meshPointsOnThisSegment
                            )
                        );

                        forAll(realCellsFoundInRange,cFIR)
                        {
                            const label realCell = realCellsFoundInRange[cFIR];

                            referredCell cellToRefer
                            (
                                mesh_,
                                Pstream::myProcNo(),
                                realCell,
                                patch.faceCentres()[faceT],
                                neighbFaceCentres[faceT],
                                patch.faceNormals()[faceT],
                                neighbFaceAreas[faceT]
                               /(mag(neighbFaceAreas[faceT]) + VSMALL)
                            );

                            referredCellsToTransfer.append(cellToRefer);

                            // add real cells found in range of processor patch
                            // to whole mesh list

                            if
                            (
                                findIndex
                                (
                                    realCellsWithinRCutMaxOfAnyReferringPatch,
                                    realCell
                                )
                             == -1
                            )
                            {
                                realCellsWithinRCutMaxOfAnyReferringPatch.append
                                (
                                    realCell
                                );
                            }
                        }
                    }

                    referredInteractionList.shrink();

                    labelList referredCellsFoundInRange
                    (
                        referredCellsInRangeOfSegment
                        (
                            referredInteractionList,
                            meshFacesOnThisSegment,
                            meshEdgesOnThisSegment,
                            meshPointsOnThisSegment
                        )
                    );

                    forAll(referredCellsFoundInRange,cFIR)
                    {
                        referredCell& existingRefCell
                            = referredInteractionList[referredCellsFoundInRange[cFIR]];

                        referredCell cellToReRefer = existingRefCell.reRefer
                        (
                            patch.faceCentres()[faceT],
                            neighbFaceCentres[faceT],
                            patch.faceNormals()[faceT],
                            neighbFaceAreas[faceT]
                           /(mag(neighbFaceAreas[faceT]) + VSMALL)
                        );

                        referredCellsToTransfer.append(cellToReRefer);
                    }
                }
            }

            referredCellsToTransfer.shrink();

            // Send these cells to the neighbouring processor.

            {
                OPstream toNeighbProc
                (
                    Pstream::blocking,
                    patch.neighbProcNo()
                );

//                 Pout << "Processor " << Pstream::myProcNo()
//                     << " sending to " << patch.neighbProcNo()
//                     << nl << referredCellsToTransfer << endl;

                toNeighbProc << referredCellsToTransfer;
            }
        }

        forAll(procPatches,pP)
        {
            const processorPolyPatch& patch =
                refCast<const processorPolyPatch>
            (
                mesh_.boundaryMesh()[procPatches[pP]]
            );


            // Receive the cells from neighbour

            List<referredCell> referredCellsFromNeighbour(patch.size());

            {
                IPstream fromNeighbProc
                (
                    Pstream::blocking,
                    patch.neighbProcNo()
                );

                fromNeighbProc >> referredCellsFromNeighbour;

//                 Pout << "Processor " << Pstream::myProcNo()
//                     << " receiving from " << patch.neighbProcNo()
//                     << nl << referredCellsFromNeighbour << endl;
            }

            // Check to see if they are duplicates, if not append
            // them to the referredInteractionList

            forAll(referredCellsFromNeighbour,rCFN)
            {
                referredCell& cellToRefer
                    = referredCellsFromNeighbour[rCFN];

                // Test all existing referred and real cells to check duplicates
                // are not being made or cells aren't being referred back onto
                // themselves

                bool addCellToRefer = true;

                // Check if cellToRefer is an existing referred cell

                forAll(referredInteractionList, rIL)
                {
                    if(cellToRefer.duplicate(referredInteractionList[rIL]))
                    {
                        addCellToRefer = false;

                        break;
                    }
                }

                // Check for cellToRefer being referred back ontop of a real cell

                if
                (
                    cellToRefer.duplicate
                    (
                        Pstream::myProcNo(),
                        mesh_.nCells()
                    )
                )
                {
                    addCellToRefer = false;
                }

                if (addCellToRefer)
                {
                    referredInteractionList.append(cellToRefer);
                }
            }
        }
    }

    if (iterationNo == 0)
    {
        // record all real cells in range of any referring patch (cyclic or
        // processor) on the first iteration when the real cells are evaluated.

        realCellsWithinRCutMaxOfAnyReferringPatch_.transfer
        (
            realCellsWithinRCutMaxOfAnyReferringPatch.shrink()
        );

        // construct {faces, edges, points}WithinRCutMaxOfAnyReferringPatch_

        DynamicList<label> realFacesWithinRCutMaxOfAnyReferringPatch;

        DynamicList<label> realEdgesWithinRCutMaxOfAnyReferringPatch;

        DynamicList<label> realPointsWithinRCutMaxOfAnyReferringPatch;

        forAll(realCellsWithinRCutMaxOfAnyReferringPatch_, rCWR)
        {
            const label realCell
            (
                realCellsWithinRCutMaxOfAnyReferringPatch_[rCWR]
            );

            const labelList& rCFaces
            (
                mesh_.cells()[realCell]
            );

            forAll(rCFaces, rCF)
            {
                const label f(rCFaces[rCF]);

                if
                (
                    findIndex(realFacesWithinRCutMaxOfAnyReferringPatch,f)
                 == -1
                )
                {
                    realFacesWithinRCutMaxOfAnyReferringPatch.append(f);
                }
            }

            const labelList& rCEdges
            (
                mesh_.cellEdges()[realCell]
            );

            forAll(rCEdges, rCE)
            {
                const label e(rCEdges[rCE]);

                if
                (
                    findIndex(realEdgesWithinRCutMaxOfAnyReferringPatch,e)
                 == -1
                )
                {
                    realEdgesWithinRCutMaxOfAnyReferringPatch.append(e);
                }
            }

            const labelList& rCPoints
            (
                mesh_.cellPoints()[realCell]
            );

            forAll(rCPoints, rCP)
            {
                const label p(rCPoints[rCP]);

                if
                (
                    findIndex(realPointsWithinRCutMaxOfAnyReferringPatch,p)
                 == -1
                )
                {
                    realPointsWithinRCutMaxOfAnyReferringPatch.append(p);
                }
            }
        }

        realFacesWithinRCutMaxOfAnyReferringPatch_.transfer
        (
            realFacesWithinRCutMaxOfAnyReferringPatch.shrink()
        );

        realEdgesWithinRCutMaxOfAnyReferringPatch_.transfer
        (
            realEdgesWithinRCutMaxOfAnyReferringPatch.shrink()
        );

        realPointsWithinRCutMaxOfAnyReferringPatch_.transfer
        (
            realPointsWithinRCutMaxOfAnyReferringPatch.shrink()
        );
    }

    iterationNo++;

    cellsReferredThisIteration =
        referredInteractionList.size() - refIntListStartSize;

//     Pout << tab << "Cells added this iteration: "
//         << cellsReferredThisIteration << endl;

    reduce(cellsReferredThisIteration, sumOp<label>());

    Info << tab << "Cells added this iteration: "
        << cellsReferredThisIteration << endl;
}

referredInteractionList.shrink();

//     Info << "referredInteractionList.size() = "
//         << referredInteractionList.size() << endl;

//     forAll(referredInteractionList,rIL)
//     {
//         Info << referredInteractionList[rIL];
//     }
