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

#include "topoChangerFvMesh.H"
#include "domainDecomposition.H"
#include "fvFieldDecomposer.H"
#include "labelIOField.H"
#include "processorMeshesReconstructor.H"
#include "fvFieldReconstructor.H"
#include "passiveProcessorPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "cloud.H"
#include "cloudDistribute.H"
#include "meshObjectBase.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::topoChangerFvMesh::loadBalance(const dictionary& decompDict)
{
    if (!Pstream::parRun())
    {
        InfoIn("bool topoChangerFvMesh::loadBalance()")
            << "Load balancing a serial run.  Skipping."
            << endl;

        return false;
    }

    // Check imbalance.  Note: add run-time selection for the imbalance
    // weights criterion

    // Create a parallel decomposition
    domainDecomposition meshDecomp
    (
        *this,
        decompDict
    );

    // Decompose the mesh based on processor data
    meshDecomp.decomposeMesh(false);

    // Analyse mesh decomposition to see how many cells are passed
    // to which processor
    labelListList migratedCells(meshDecomp.nProcs());

    // Fill local list with number of local cells to be sent to each processor
    // First index = my processor.  Second index = processor to send to
    labelList& curMigratedCells = migratedCells[Pstream::myProcNo()];

    curMigratedCells.setSize(Pstream::nProcs(), 0);

    const labelList& cellToProc = meshDecomp.cellToProc();

    // // Write as volScalarField for post-processing
    // volScalarField cellDist
    // (
    //     IOobject
    //     (
    //         "cellDist",
    //         time().timeName(),
    //         *this,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     *this,
    //     dimensionedScalar("cellDist", dimless, 0),
    //     zeroGradientFvPatchScalarField::typeName
    // );

    // forAll(cellToProc, celli)
    // {
    //     cellDist[celli] = cellToProc[celli];
    // }

    // cellDist.write();

    // Count migrated cells
    forAll (cellToProc, cellI)
    {
        curMigratedCells[cellToProc[cellI]]++;
    }

    // Now that each processor has filled in its own part, combine the data
    Pstream::gatherList(migratedCells);
    Pstream::scatterList(migratedCells);

    if (debug)
    {
        Info<< "Migrated cells per processor: " << migratedCells << endl;
    }

    // Reading through second index now tells how many cells will arrive
    // from which processor
    forAll (migratedCells, i)
    {
        bool allZero = true;
        forAll (migratedCells, j)
        {
            if (migratedCells[j][i] > 0)
            {
                allZero = false;
                break;
            }
        }

        if (allZero)
        {
            FatalErrorIn("bool topoChangerFvMesh::loadBalance()")
                << "Reconstructed mesh must have at least one cell "
                    "on each processor"
                << abort(FatalError);
        }
    }

    // Find out which processor faces will become local internal faces
    // by comparing decomposition index from the other side


    // Split the mesh and fields over processors and send

    // Prepare receiving side

    // Create the reconstructor

    // HR 21.12.18 : Use empty domainname to avoid auto-created of
    // fvSchemes/fvSolution
    processorMeshesReconstructor meshRecon("");

    PtrList<fvMesh>& procMeshes = meshRecon.meshes();
    procMeshes.setSize(meshDecomp.nProcs());
    meshRecon.globalPointIndex().setSize(meshDecomp.nProcs());

    // Collect local fields for decomposition
    clearOut();

    // Collect fields for load balancing
    HashTable<const volScalarField*> volScalarFields =
        thisDb().lookupClass<volScalarField>();

    HashTable<const volVectorField*>volVectorFields =
        thisDb().lookupClass<volVectorField>();

    HashTable<const surfaceScalarField*> surfaceScalarFields =
        thisDb().lookupClass<surfaceScalarField>();

    HashTable<const surfaceVectorField*> surfaceVectorFields =
        thisDb().lookupClass<surfaceVectorField>();

    // Particles
    HashTable<const cloud*> clouds = thisDb().lookupClass<cloud>();

    // Distribute cell and point level for AMR + DLB runs. VV, 18/May/2018

    // Check cellLevel and pointLevel
    // Note: possible sync problem if cellLevel or pointLevel is not present
    // on all processors in comms.  However, parallel reduce check
    // is not allowed
    // HJ, 22/Oct/2018
    bool cellLevelFound = this->foundObject<labelIOField>("cellLevel");
    bool pointLevelFound = this->foundObject<labelIOField>("pointLevel");

    //HJ, HERE: remove the fields that should not be load balanced

    // Prepare fields for reconstruction
    // First index: field index from volScalarFields
    // Second index: processor entry for the field
    List<PtrList<volScalarField> > receivedVolScalarFields
    (
        volScalarFields.size()
    );

    List<PtrList<volVectorField> > receivedVolVectorFields
    (
        volVectorFields.size()
    );

    List<PtrList<surfaceScalarField> > receivedSurfaceScalarFields
    (
        surfaceScalarFields.size()
    );

    List<PtrList<surfaceVectorField> > receivedSurfaceVectorFields
    (
        surfaceVectorFields.size()
    );

    forAll (receivedVolScalarFields, fieldI)
    {
        receivedVolScalarFields[fieldI].setSize(Pstream::nProcs());
    }

    forAll (receivedVolVectorFields, fieldI)
    {
        receivedVolVectorFields[fieldI].setSize(Pstream::nProcs());
    }

    forAll (receivedSurfaceScalarFields, fieldI)
    {
        receivedSurfaceScalarFields[fieldI].setSize(Pstream::nProcs());
    }

    forAll (receivedSurfaceVectorFields, fieldI)
    {
        receivedSurfaceVectorFields[fieldI].setSize(Pstream::nProcs());
    }

    // Cell and point level
    // Note: ordinary lists instead of IO lists
    PtrList<labelList> receivedCellLevel(Pstream::nProcs());
    PtrList<labelList> receivedPointLevel(Pstream::nProcs());

    // Clouds
    PtrList<cloudDistribute> cloudDistributes(clouds.size());
    {
        label cloudI = 0;
        forAllConstIter(HashTable<const cloud*>, clouds, iter)
        {
            cloud& c = const_cast<cloud&>(*iter());
            cloudDistributes.set
            (
                cloudI++,
                c.cloudDist
                (
                    meshDecomp.cellToProc(),
                    meshDecomp.procCellAddressing(),
                    meshDecomp.procFaceAddressing()
                )
            );
        }
    }

    // HR 14.12.18: Create dummy database pointing into the non-parallel case
    // directory and disable the runTimeModifiable switch.  dummyTime is used
    // for decomposed, received and reconstructed fvMeshes (ie. before its data
    // is moved into *this).
    //
    // The pointing into the non-parallel case directory is somewhat a hack to
    // prevent auto-creation of fvSchemes and fvSolution (they exist in the root
    // case). This is potentially dangerous if we write anything, but fvMeshes
    // derived from this database are NO_WRITE so we should be fine. Making if
    // non-runTimeModifiable prevents registration of fvSolution with the
    // fileMonitor (addWatch). Not all processors will necessarily receive a
    // mesh and watches will then cause dead-locks!
    Time dummyTime
    (
        time().rootPath(),
        time().globalCaseName(),
        "system",
        "constant",
        false
    );

    const_cast<Switch&>(dummyTime.runTimeModifiable()) = false;

    for (label procI = 0; procI < meshDecomp.nProcs(); procI++)
    {
        // Check if there is a mesh to send
        if (migratedCells[Pstream::myProcNo()][procI] > 0)
        {
            // Send mesh and field to processor procI

            // Get mesh from decomposition
            autoPtr<fvMesh> procMeshPtr = meshDecomp.processorMesh
            (
                procI,
                dummyTime,
                "",     // HR 21.12.18 : Use empty domainname to avoid
                        // auto-created offvSchemes/fvSolution
                true    // Create passive processor patches
            );
            fvMesh& procMesh = procMeshPtr();

            // Create a field decomposer
            fvFieldDecomposer fieldDecomposer
            (
                *this,
                procMesh,
                meshDecomp.procFaceAddressing()[procI],
                meshDecomp.procCellAddressing()[procI],
                meshDecomp.procBoundaryAddressing()[procI]
            );

            if (procI != Pstream::myProcNo())
            {
                if (debug)
                {
                    Pout<< "Send mesh and fields to processor " << procI
                        << endl;
                }

                OPstream toProc
                (
                    Pstream::blocking,
                    procI
                );

                // Send the mesh and fields to target processor
                toProc << procMesh << nl;
                toProc << meshDecomp.globalPointIndex(procI) << nl;

                // Send fields
                sendFields(volScalarFields, fieldDecomposer, toProc);
                sendFields(volVectorFields, fieldDecomposer, toProc);
                sendFields(surfaceScalarFields, fieldDecomposer, toProc);
                sendFields(surfaceVectorFields, fieldDecomposer, toProc);

                // Send cell level with procCellAddressing
                if (cellLevelFound)
                {
                    const labelIOField& cellLevel =
                        this->lookupObject<labelIOField>("cellLevel");

                    toProc <<
                        labelList
                        (
                            cellLevel,
                            meshDecomp.procCellAddressing()[procI]
                        )
                        << nl;
                }

                if (pointLevelFound)
                {
                    const labelIOField& pointLevel =
                        this->lookupObject<labelIOField>("pointLevel");

                    // Send point level with procPointAddressing
                    toProc <<
                        labelList
                        (
                            pointLevel,
                            meshDecomp.procPointAddressing()[procI]
                        )
                        << nl;
                }

                // Send clouds
                forAll(cloudDistributes, cloudI)
                {
                    cloudDistributes[cloudI].send(toProc, procI);
                }
            }
            else
            {
                // My processor.  Grab mesh and fields
                // Note: procI == Pstream::myProcNo()
                // HJ, 26/Apr/2018
                // Pout<< "Set local mesh and fields" << endl;
                // Set local mesh piece
                procMeshes.set
                (
                    procI,
                    procMeshPtr
                );

                meshRecon.globalPointIndex()[procI] =
                    meshDecomp.globalPointIndex(procI);

                // Set local fields
                // Note: first index is field index and second index is procI
                insertFields
                (
                    volScalarFields,
                    fieldDecomposer,
                    receivedVolScalarFields
                );

                insertFields
                (
                    volVectorFields,
                    fieldDecomposer,
                    receivedVolVectorFields
                );

                insertFields
                (
                    surfaceScalarFields,
                    fieldDecomposer,
                    receivedSurfaceScalarFields
                );

                insertFields
                (
                    surfaceVectorFields,
                    fieldDecomposer,
                    receivedSurfaceVectorFields
                );

                // Insert cell level
                if (cellLevelFound)
                {
                    const labelIOField& cellLevel =
                        this->lookupObject<labelIOField>("cellLevel");

                    // Insert cell level
                    receivedCellLevel.set
                    (
                        Pstream::myProcNo(),
                        new labelList
                        (
                            cellLevel,
                            meshDecomp.procCellAddressing()
                                [Pstream::myProcNo()]
                        )
                    );
                }

                // Insert point level
                if (pointLevelFound)
                {
                    const labelIOField& pointLevel =
                        this->lookupObject<labelIOField>("pointLevel");

                    // Insert point level
                    receivedPointLevel.set
                    (
                        Pstream::myProcNo(),
                        new labelList
                        (
                            pointLevel,
                            meshDecomp.procPointAddressing()
                                [Pstream::myProcNo()]
                        )
                    );
                }

                // HR, 18.11.2018. Distribution of clouds is trivial and is
                // treated in Cloud<ParticleType>::split, which is called in
                // the constructor of CloudDistribute.
            }
        }
    }

    // Collect pieces of mesh and fields from other processors
    for (label procI = 0; procI < meshDecomp.nProcs(); procI++)
    {
        if (procI != Pstream::myProcNo())
        {
            // Check if there is a mesh to send
            if (migratedCells[procI][Pstream::myProcNo()] > 0)
            {
                if (debug)
                {
                    Pout<< "Receive mesh and fields from " << procI << endl;
                }

                // Note: communication can be optimised.  HJ, 27/Feb/2018
                IPstream fromProc
                (
                    Pstream::blocking,
                    procI
                );

                // Receive the mesh
                procMeshes.set
                (
                    procI,
                    new fvMesh
                    (
                        IOobject
                        (
                            "", // "processorPart" + Foam::name(procI),
                            dummyTime.timeName(),
                            dummyTime,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        fromProc,
                        false             // Do not sync
                    )
                );

                // Receive the global points addr
                labelList gppi(fromProc);
                meshRecon.globalPointIndex()[procI] = gppi;

                // Receive fields
                // Note: first index is field index and second index is procI
                receiveFields
                (
                    procI,
                    receivedVolScalarFields,
                    procMeshes[procI],
                    fromProc
                );

                receiveFields
                (
                    procI,
                    receivedVolVectorFields,
                    procMeshes[procI],
                    fromProc
                );

                receiveFields
                (
                    procI,
                    receivedSurfaceScalarFields,
                    procMeshes[procI],
                    fromProc
                );

                receiveFields
                (
                    procI,
                    receivedSurfaceVectorFields,
                    procMeshes[procI],
                    fromProc
                );

                // Receive cell level
                if (cellLevelFound)
                {
                    receivedCellLevel.set
                    (
                        procI,
                        new labelList(fromProc)
                    );
                }

                // Receive point level
                if (pointLevelFound)
                {
                    receivedPointLevel.set
                    (
                        procI,
                        new labelList(fromProc)
                    );
                }

                // Receive clouds
                forAll(cloudDistributes, cloudI)
                {
                    cloudDistributes[cloudI].receive(fromProc, procI);
                }
            }
        }
    }

    if (debug)
    {
        forAll (procMeshes, procI)
        {
            if (procMeshes.set(procI))
            {
                Pout<< "procMesh " << procI
                    << " points " << procMeshes[procI].nPoints()
                    << " faces: " << procMeshes[procI].nFaces()
                    << " internal: " << procMeshes[procI].nInternalFaces()
                    << " cells: " << procMeshes[procI].nCells()
                    << " patches: " << procMeshes[procI].boundary().size()
                    << endl;
            }
        }
    }

    // Create the reconstructed mesh
    autoPtr<fvMesh> reconstructedMeshPtr =
        meshRecon.reconstructMesh(dummyTime);
    fvMesh& reconMesh = reconstructedMeshPtr();

    if (debug)
    {
        Pout<< "Reconstructed mesh stats: "
            << " nCells: " << reconMesh.nCells()
            << " nFaces: " << reconMesh.nFaces()
            << " nIntFaces: " << reconMesh.nInternalFaces()
            << " polyPatches: "
            << reconMesh.boundaryMesh().size()
            << " patches: "
            << reconMesh.boundary().size()
            << endl;
    }

    // Apply changes to the local mesh:
    // - refactor the boundary to match new patches.  Note: processor
    //   patch types may be added or removed
    // - reset all primitives

    // Use non-const access to boundaryMesh
    polyBoundaryMesh& bMesh = const_cast<polyBoundaryMesh&>(boundaryMesh());

    const polyBoundaryMesh& reconBMesh = reconMesh.boundaryMesh();

    // Resize the existing boundary to match in the number of patches
    if (bMesh.size() > reconBMesh.size())
    {
        // Decreasing in size: check patches that are due to disappear

        // Check before resizing: only processor patches may be deleted
        bool okToDelete = true;

        for (label patchI = reconBMesh.size(); patchI < bMesh.size(); patchI++)
        {
            if (!isA<processorPolyPatch>(bMesh[patchI]))
            {
                okToDelete = false;
                break;
            }
        }

        if (!okToDelete)
        {
            FatalErrorIn("bool topoChangerFvMesh::loadBalance()")
                << "Boundary resizing error: deleting non-processor patches "
                << bMesh.size() << " to " << reconBMesh.size() << nl
                << "old patch types: " << bMesh.types() << nl
                << "new patch types: " << reconBMesh.types() << nl
                << "This is not allowed."
                << abort(FatalError);
        }

        // Resize the boundary.  Patches hanging off the end of the list
        // will be deleted by PtrList
        bMesh.setSize(reconBMesh.size());
    }
    else if (bMesh.size() < reconBMesh.size())
    {
        // Increasing in size: check patches that are due to disappear
        bMesh.setSize(reconBMesh.size());
    }

    // Delete processor patches from old boundary
    boolList patchFlag(bMesh.size(), false);

    forAll (bMesh, patchI)
    {
        // If the patch is not set, the slot is available for re-use
        if (!bMesh.set(patchI))
        {
            patchFlag[patchI] = true;
        }
        else if (isA<processorPolyPatch>(bMesh[patchI]))
        {
            // If the patch is set and contains an "old" processor patch
            // the old patch will be deleted and replaced with a new one
            patchFlag[patchI] = true;
        }
    }

    // Store the resetFvPatchFlag to indicate patches that will be rebuilt
    // patchFlag will be used in further signalling
    boolList resetFvPatchFlag = patchFlag;

    // Check alignment and types of old and new boundary
    // Insert new patches processor patches
    // Collect patch sizes and starts

    labelList reconPatchSizes(reconBMesh.size(), 0);

    forAll (bMesh, patchI)
    {
        // If the patch is preserved, types need to match
        if (!patchFlag[patchI])
        {
            const polyPatch& bPatch = bMesh[patchI];
            const polyPatch& reconBPatch = reconBMesh[patchI];

            // Check if the patch is matching
            if
            (
                bPatch.type() != reconBPatch.type()
             || bPatch.name() != reconBPatch.name()
            )
            {
                FatalErrorIn("bool topoChangerFvMesh::loadBalance()")
                    << "Patch names or types not matching for index "
                    << patchI << nl
                    << "bMesh name: " << bPatch.name()
                    << " type: " << bPatch.type() << nl
                    << "reconBMesh name: " << reconBPatch.name()
                    << "  type: " << reconBPatch.type() << nl
                    << abort(FatalError);
            }

            // Record patch size
            reconPatchSizes[patchI] = reconBMesh[patchI].size();
        }
    }

    // Insert new processor patches
    // Note:
    // The code is set up so that any intermediate empty boundary patch
    // slots are re-used.  It is expected that the new patches will be
    // added at the end, so this is safety.
    // As a result, patch sizes and starts need to be re-assembled
    // rather than used directly from reconBMesh
    // HJ, 16/Apr/2018
    forAll (reconBMesh, reconPatchI)
    {
        if (isA<processorPolyPatch>(reconBMesh[reconPatchI]))
        {
            // Find a slot to insert into
            label nextPatchSlot = 0;

            for
            (
                nextPatchSlot = 0;
                nextPatchSlot < patchFlag.size();
                nextPatchSlot++
            )
            {
                if (patchFlag[nextPatchSlot])
                {
                    // Found a patch slot
                    break;
                }
            }

            if (nextPatchSlot >= patchFlag.size())
            {
                FatalErrorIn("bool topoChangerFvMesh::loadBalance()")
                    << "Cannot find available patch slot: " << nextPatchSlot
                    << abort(FatalError);
            }

            // Insert a patch
            const processorPolyPatch& ppPatch =
                refCast<const processorPolyPatch>
                (
                    reconBMesh[reconPatchI]
                );

            // Note:
            // The new mesh has not been reset.  Therefore, sizes and starts
            // cannot be set properly until the mesh is reset
            // This is done in the resetPrimitives function
            // HJ, 16/Apr/2018
            bMesh.set
            (
                nextPatchSlot,
                new processorPolyPatch
                (
                    // ppPatch.name(),
                    "procBoundary" + Foam::name(ppPatch.myProcNo()) + "to"
                    + Foam::name(ppPatch.neighbProcNo()),
                    0,                  // dummy size
                    nInternalFaces(),   // dummy start
                    // ppPatch.size(),
                    // ppPatch.start(),
                    nextPatchSlot,
                    bMesh,
                    ppPatch.myProcNo(),
                    ppPatch.neighbProcNo()
                )
            );

            // Mark patch slot as used
            patchFlag[nextPatchSlot] = false;

            // Collect size
            reconPatchSizes[nextPatchSlot] = reconBMesh[reconPatchI].size();
        }
    }

    // Patch sizes are assembled.  Collect patch starts for the new mesh
    labelList reconPatchStarts(reconBMesh.size(), 0);

    if (!reconPatchStarts.empty())
    {
        reconPatchStarts[0] = reconMesh.nInternalFaces();

        for (label i = 1; i < reconPatchStarts.size(); i++)
        {
            reconPatchStarts[i] =
                reconPatchStarts[i - 1]
              + reconPatchSizes[i - 1];
        }
    }

    // Create mesh map.  Note: map is dummy and it is used only for resizing
    // HJ, 16/May/2018
    const polyBoundaryMesh& patches = boundaryMesh();
    labelList patchStarts(patches.size());
    labelList oldPatchNMeshPoints(patches.size());
    labelListList patchPointMap(patches.size());

    forAll(patches, patchI)
    {
        patchStarts[patchI] = patches[patchI].start();
        oldPatchNMeshPoints[patchI] = patches[patchI].nPoints();
        patchPointMap[patchI].setSize(patches[patchI].nPoints(), -1);
    }

    mapPolyMesh meshMap
    (
        *this,                 // const polyMesh& mesh,
        nPoints(),             // nOldPoints,
        nFaces(),              // nOldFaces,
        nCells(),              // nOldCells,

        labelList(reconMesh.nPoints(), -1),   // pointMap,
        List<objectMap>(0),         // pointsFromPoints,

        labelList(reconMesh.nFaces(), -1),   // faceMap,
        List<objectMap>(0),         // facesFromPoints,
        List<objectMap>(0),         // facesFromEdges,
        List<objectMap>(0),         // facesFromFaces,

        labelList(),                // cellMap,
        List<objectMap>(0),         // cellsFromPoints,
        List<objectMap>(0),         // cellsFromEdges,
        List<objectMap>(0),         // cellsFromFaces,
        List<objectMap>(0),         // cellsFromCells,

        labelList(nPoints(), -1),   // reversePointMap,
        labelList(nFaces(), -1),    // reverseFaceMap,
        labelList(nCells(), -1),    // reverseCellMap,

        labelHashSet(0),            // flipFaceFlux,

        labelListList(0),           // patchPointMap,
        labelListList(0),           // pointZoneMap,
        labelListList(0),           // faceZonePointMap,
        labelListList(0),           // faceZoneFaceMap,
        labelListList(0),           // cellZoneMap,

        resetFvPatchFlag,           // resetPatchFlag

        pointField(0),              // preMotionPoints,
        patchStarts,                // oldPatchStarts,
        oldPatchNMeshPoints         // oldPatchNMeshPoints
    );

    // Reset fvMesh and patches
    resetFvPrimitives
    (
        xferCopy(reconMesh.allPoints()),
        xferCopy(reconMesh.allFaces()),
        xferCopy(reconMesh.faceOwner()),
        xferCopy(reconMesh.faceNeighbour()),
        reconPatchSizes,
        reconPatchStarts,
        resetFvPatchFlag,
        true                       // Valid boundary
    );

    // Create field reconstructor
    fvFieldReconstructor fieldReconstructor
    (
        *this,
        procMeshes,
        meshRecon.faceProcAddressing(),
        meshRecon.cellProcAddressing(),
        meshRecon.boundaryProcAddressing()
    );

    // Rebuild fields

    rebuildFields
    (
        volScalarFields,
        fieldReconstructor,
        receivedVolScalarFields,
        meshMap
    );

    rebuildFields
    (
        volVectorFields,
        fieldReconstructor,
        receivedVolVectorFields,
        meshMap
    );

    rebuildFields
    (
        surfaceScalarFields,
        fieldReconstructor,
        receivedSurfaceScalarFields,
        meshMap
    );

    rebuildFields
    (
        surfaceVectorFields,
        fieldReconstructor,
        receivedSurfaceVectorFields,
        meshMap
    );

    // Rebuild cell level field from components
    if (cellLevelFound)
    {
        // Get (non-const) reference to cellLevel
        labelIOField& cellLevel = const_cast<labelIOField&>
            (this->lookupObject<labelIOField>("cellLevel"));

        if (cellLevel.size() != this->nCells())
        {
            cellLevel.setSize(this->nCells());
        }

        forAll (receivedCellLevel, procI)
        {
            if (receivedCellLevel.set(procI))
            {
                cellLevel.rmap
                (
                    receivedCellLevel[procI],
                    meshRecon.cellProcAddressing()[procI]
                );
            }
        }
    }

    // Rebuild point level field from components
    if (pointLevelFound)
    {
        // Get (non-const) reference to pointLevel
        labelIOField& pointLevel = const_cast<labelIOField&>
            (this->lookupObject<labelIOField>("pointLevel"));

        if (pointLevel.size() != this->nPoints())
        {
            pointLevel.setSize(this->nPoints());
        }

        forAll (receivedPointLevel, procI)
        {
            if (receivedPointLevel.set(procI))
            {
                pointLevel.rmap
                (
                    receivedPointLevel[procI],
                    meshRecon.pointProcAddressing()[procI]
                );
            }
        }
    }

    // Rebuild clouds
    forAll(cloudDistributes, cloudI)
    {
        cloudDistributes[cloudI].rebuild
        (
            meshRecon.cellProcAddressing(),
            meshRecon.faceProcAddressing()
        );
    }

    // HR 13.12.18: Update the mesh objects
    meshObjectBase::allUpdateTopology<polyMesh>(*this, meshMap);

    if (debug)
    {
        Info<< "Checking reconstructed mesh after load balancing..." << endl;
        checkMesh(true);
    }

    return true;
}


// ************************************************************************* //
