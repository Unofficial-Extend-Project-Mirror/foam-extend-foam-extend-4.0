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

#include "processorMeshesReconstructor.H"
#include "processorPolyPatch.H"
#include "passiveProcessorPolyPatch.H"
#include "SortableList.H"
#include "sharedPoints.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::processorMeshesReconstructor::firstValidMesh() const
{
    forAll (meshes_, procI)
    {
        if (meshes_.set(procI))
        {
            return procI;
        }
    }

    FatalErrorIn("label processorMeshesReconstructor::firstValidMesh() const")
        << "Cannot find a valid mesh in reconstruction set"
        << abort(FatalError);

    return 0;
}


bool Foam::processorMeshesReconstructor::readMapping()
{
    // Check for mapping
    Info<< "Check for mesh mapping data for instance.  ";

    bool readOk = true;

    forAll (meshes_, procI)
    {
        if (meshes_.set(procI))
        {
            const fvMesh& procMesh = meshes_[procI];

            IOobject pointProcAddressingHeader
            (
                "pointProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::MUST_READ
            );

            IOobject faceProcAddressingHeader
            (
                "faceProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            );

            IOobject cellProcAddressingHeader
            (
                "cellProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            );

            IOobject boundaryProcAddressingHeader
            (
                "boundaryProcAddressing",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            );

            if
            (
                !pointProcAddressingHeader.headerOk()
             || !faceProcAddressingHeader.headerOk()
             || !cellProcAddressingHeader.headerOk()
             || !boundaryProcAddressingHeader.headerOk()
            )
            {
                readOk = false;
                break;
            }
        }
    }

    // All processors are fine: read mapping data
    if (readOk)
    {
        Info<< "Mapping data present.  Reading." << endl;

        // Size the mapping arrays
        pointProcAddressing_.setSize(meshes_.size());
        faceProcAddressing_.setSize(meshes_.size());
        cellProcAddressing_.setSize(meshes_.size());
        boundaryProcAddressing_.setSize(meshes_.size());

        forAll (meshes_, procI)
        {
            if (meshes_.set(procI))
            {
                const fvMesh& procMesh = meshes_[procI];

                pointProcAddressing_.set
                (
                    procI,
                    new labelIOList
                    (
                        IOobject
                        (
                            "pointProcAddressing",
                            procMesh.facesInstance(),
                            procMesh.meshSubDir,
                            procMesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        )
                    )
                );

                faceProcAddressing_.set
                (
                    procI,
                    new labelIOList
                    (
                        IOobject
                        (
                            "faceProcAddressing",
                            procMesh.facesInstance(),
                            procMesh.meshSubDir,
                            procMesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        )
                    )
                );

                cellProcAddressing_.set
                (
                    procI,
                    new labelIOList
                    (
                        IOobject
                        (
                            "cellProcAddressing",
                            procMesh.facesInstance(),
                            procMesh.meshSubDir,
                            procMesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        )
                    )
                );

                boundaryProcAddressing_.set
                (
                    procI,
                    new labelIOList
                    (
                        IOobject
                        (
                            "boundaryProcAddressing",
                            procMesh.facesInstance(),
                            procMesh.meshSubDir,
                            procMesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        )
                    )
                );
            }
        }

        Info<< "Addressing from files: " << endl;
        forAll (meshes_, procI)
        {
            if (meshes_.set(procI))
            {
                Info<< "Proc " << procI
                    << " point addr: " << pointProcAddressing_[procI].size()
                    << " face addr: " << faceProcAddressing_[procI].size()
                    << " cell addr: " << cellProcAddressing_[procI].size()
                    << " boundary addr: "
                    << boundaryProcAddressing_[procI].size()
                    << endl;
            }
        }
    }
    else
    {
        Info<< "No mapping data available." << endl;
    }

    return readOk;
}


void Foam::processorMeshesReconstructor::writeAddressing()
{
    forAll (pointProcAddressing_, procI)
    {
        if (meshes_.set(procI))
        {
            pointProcAddressing_[procI].write();
            faceProcAddressing_[procI].write();
            cellProcAddressing_[procI].write();
            boundaryProcAddressing_[procI].write();
        }
    }
}


const Foam::processorPolyPatch&
Foam::processorMeshesReconstructor::neighbourProcPatch
(
    const processorPolyPatch& procPatch
) const
{
    const label masterProcID = procPatch.neighbProcNo();

    if (!meshes_.set(masterProcID))
    {
        FatalErrorIn
        (
            "const processorPolyPatch&\n"
            "processorMeshesReconstructor::neighbourProcPatch\n"
            "(\n"
            "    const processorPolyPatch& procPatch\n"
            ") const"
        )   << "Cannot find processor patch pair ("
            << procPatch.myProcNo() << " "
            << procPatch.neighbProcNo() << ") for merging"
            << abort(FatalError);
    }

    const polyMesh& masterMesh = meshes_[masterProcID];

    bool found = false;

    // Find the processor patch that corresponds to current patch
    const polyBoundaryMesh& masterPatches = masterMesh.boundaryMesh();

    forAll (masterPatches, masterPatchI)
    {
        if
        (
            isA<processorPolyPatch>
            (
                masterPatches[masterPatchI]
            )
        )
        {
            const processorPolyPatch& masterProcPatch =
                refCast<const processorPolyPatch>
                (
                    masterPatches[masterPatchI]
                );

            // Check neighbour processor index
            if (masterProcPatch.neighbProcNo() == procPatch.myProcNo())
            {
                // Found matching patch.  Check sizes
                if (masterProcPatch.size() != procPatch.size())
                {
                    FatalErrorIn
                    (
                        "const processorPolyPatch&\n"
                        "processorMeshesReconstructor::neighbourProcPatch\n"
                        "(\n"
                        "    const processorPolyPatch& procPatch\n"
                        ") const"
                    )   << "Processor patch pair ("
                        << procPatch.myProcNo() << " "
                        << procPatch.neighbProcNo() << ") sizes do not match: "
                        << masterProcPatch.size() << " and " << procPatch.size()
                        << abort(FatalError);
                }

                return masterProcPatch;
            }
        }
    }

    if (!found)
    {
        FatalErrorIn
        (
            "const processorPolyPatch&\n"
            "processorMeshesReconstructor::neighbourProcPatch\n"
            "(\n"
            "    const processorPolyPatch& procPatch\n"
            ") const"
        )   << "Cannot find processor patch pair ("
            << procPatch.myProcNo() << " "
            << procPatch.neighbProcNo() << ") for merging"
            << abort(FatalError);
    }

    // Dummy return
    return procPatch;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::processorMeshesReconstructor::reconstructPoints(fvMesh& mesh) const
{
    // Create the new points
    vectorField newPoints(mesh.nPoints());

    forAll (meshes_, procI)
    {
        if (meshes_.set(procI))
        {
            // Reconstruct only live points.  HJ, 7/Mar/2011
            const vectorField& procPoints = meshes_[procI].points();

            // Set the cell values in the reconstructed field

            const labelList& pointProcAddressingI =
                pointProcAddressing()[procI];

            if (pointProcAddressingI.size() != procPoints.size())
            {
                FatalErrorIn("processorMeshes")
                    << "problem :"
                    << " pointProcAddressingI: " << pointProcAddressingI.size()
                    << " procPoints:" << procPoints.size()
                    << abort(FatalError);
            }

            // Only live points carry reconstruction data.  Reconsider
            // HJ, 6/Sep/2009
            for (label pointI = 0; pointI < meshes_[procI].nPoints(); pointI++)
            {
                newPoints[pointProcAddressingI[pointI]] = procPoints[pointI];
            }
        }
    }

    mesh.movePoints(newPoints);
    mesh.write();
}


Foam::autoPtr<Foam::fvMesh>
Foam::processorMeshesReconstructor::reconstructMesh(const Time& db)
{
    // Note:
    // In load balancing, there will exist a double set of processor patches
    // One, created by moving cells adjacent to "old" processor boundaries
    // (processorPolyPatch) and another, created by splitting up previously
    // internal mesh faces into "new" processor patches
    // In order to reconstruct the mesh correctly, the two sets of processor
    // boundaries are kept separately
    // In reconstruction, the order of processor patch faces needs to be
    // preserved.  This is achieved by
    // - first adding original processor patch faces
    // - adding passiveProcessor patch faces in the processor order

    // Check for read

    if (readMapping())
    {
        // Mapping data present and read.  Reconstructed mesh may be
        // present as well
        bool readMesh = false;

        {
            IOobject points
            (
                "points",
                db.timeName(),
                meshes_[firstValidMesh()].meshSubDir,
                db,
                IOobject::MUST_READ
            );

            if (points.headerOk())
            {
                readMesh = true;
            }
        }

        if (readMesh)
        {
            Info<< "Global mesh present for time " << db.timeName()
                << ".  Reading mesh." << endl;

            autoPtr<fvMesh> globalMeshPtr
            (
                new fvMesh
                (
                    IOobject
                    (
                        meshName_,
                        db.timeName(),
                        db,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                )
            );

            return globalMeshPtr;
        }
    }

    // Get first boundary mesh size
    const label firstBndryMeshSize =
        meshes_[firstValidMesh()].boundaryMesh().size();

    // Prepare patch reconstruction. Note: default key type is word in HashTable
    HashTable<label> patchNameLookup(firstBndryMeshSize);
    DynamicList<word> patchTypeLookup(firstBndryMeshSize);

    label nReconPatches = 0;

    forAll (meshes_, procI)
    {
        if (meshes_.set(procI))
        {
            const polyBoundaryMesh& procPatches = meshes_[procI].boundaryMesh();

            forAll (procPatches, patchI)
            {
                if (!isA<processorPolyPatch>(procPatches[patchI]))
                {
                    // Regular patch, get the patch and the name
                    const polyPatch& patch = procPatches[patchI];
                    const word& patchName = patch.name();

                    // Insert the pair into the hash table
                    if (patchNameLookup.insert(patchName, nReconPatches))
                    {
                        // This is the first time we're inserting this patch,
                        // add its type into the list and increment the counter
                        patchTypeLookup.append(patch.type());
                        ++nReconPatches;
                    }
                    // else already present in the table
                }
            }
        }
    }

    // Fill in patch names
    wordList reconPatchNames(patchNameLookup.size());

    // Loop through all the patch names and collect names and types
    forAllConstIter(HashTable<label>, patchNameLookup, iter)
    {
        reconPatchNames[iter()] = iter.key();
    }

    // Transfer the contents of the patchTypeLookup dynamic list into the
    // ordinary list. Note: the ordering remains the same
    const wordList reconPatchTypes(patchTypeLookup.xfer());

    // Prepare point, face and patch reconstruction
    label nReconPoints = 0;
    label nReconFaces = 0;
    label nReconCells = 0;
    labelList reconPatchSizes(reconPatchTypes.size(), 0);

    forAll (meshes_, procI)
    {
        if (meshes_.set(procI))
        {
            // Get current mesh
            const fvMesh& procMesh = meshes_[procI];

            // Count total number of points and faces
            nReconPoints += procMesh.nPoints();
            nReconFaces += procMesh.allFaces().size();
            nReconCells += procMesh.nCells();

            const polyBoundaryMesh& procPatches = procMesh.boundaryMesh();

            forAll (procPatches, patchI)
            {
                if (!isA<processorPolyPatch>(procPatches[patchI]))
                {
                    // Get current patch
                    const polyPatch& patch = procPatches[patchI];

                    // Find processor patch index in reconstructed boundary
                    const label pnIndex = patchNameLookup.find(patch.name())();

                    // Check patch name and type
                    if
                    (
                        patch.name() != reconPatchNames[pnIndex]
                     || patch.type() != reconPatchTypes[pnIndex]
                    )
                    {
                        FatalErrorIn
                        (
                            "autoPtr<fvMesh> processorMeshesReconstructor::"
                            "reconstructMesh(const Time& db)"
                        )   << "Patch name and type does not match "
                            << "across processors for patch "
                            << patch.name() << " type: " << patch.type()
                            << abort(FatalError);
                    }

                    // Record number of faces in patch
                    reconPatchSizes[pnIndex] += patch.size();
                }
            }
        }
    }

    // Note: for easier debugging, set mapping, owner and neighbour to -1
    pointField reconPoints(nReconPoints);
    faceList reconFaces(nReconFaces);
    labelList cellOffset(meshes_.size(), 0);
    labelList reconOwner(nReconFaces, -1);
    labelList reconNeighbour(nReconFaces, -1);
    faceListList reconPatchFaces(reconPatchTypes.size());
    labelListList reconPatchOwner(reconPatchTypes.size());

    forAll (reconPatchFaces, patchI)
    {
        // Set size of reconstructed patch faces
        reconPatchFaces[patchI].setSize(reconPatchSizes[patchI]);

        // Set size of reconstructed patch face owners with default value of -1
        reconPatchOwner[patchI].setSize(reconPatchSizes[patchI], -1);
    }

    // Size the mapping arrays
    pointProcAddressing_.setSize(meshes_.size());
    faceProcAddressing_.setSize(meshes_.size());
    cellProcAddressing_.setSize(meshes_.size());
    boundaryProcAddressing_.setSize(meshes_.size());

    // Allocate addressing arrays on all meshes
    forAll (meshes_, procI)
    {
        if (meshes_.set(procI))
        {
            const fvMesh& procMesh = meshes_[procI];

            pointProcAddressing_.set
            (
                procI,
                new labelIOList
                (
                    IOobject
                    (
                        "pointProcAddressing",
                        procMesh.facesInstance(),
                        procMesh.meshSubDir,
                        procMesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    labelList(procMesh.nPoints(), -1)
                )
            );

            faceProcAddressing_.set
            (
                procI,
                new labelIOList
                (
                    IOobject
                    (
                        "faceProcAddressing",
                        procMesh.facesInstance(),
                        procMesh.meshSubDir,
                        procMesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    labelList(procMesh.allFaces().size(), -1)
                )
            );

            cellProcAddressing_.set
            (
                procI,
                new labelIOList
                (
                    IOobject
                    (
                        "cellProcAddressing",
                        procMesh.facesInstance(),
                        procMesh.meshSubDir,
                        procMesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    labelList(procMesh.nCells(), -1)
                )
            );

            boundaryProcAddressing_.set
            (
                procI,
                new labelIOList
                (
                    IOobject
                    (
                        "boundaryProcAddressing",
                        procMesh.facesInstance(),
                        procMesh.meshSubDir,
                        procMesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    labelList(procMesh.boundaryMesh().size(), -1)
                )
            );
        }
    }

    // Reset the counters
    nReconPoints = 0;
    nReconFaces = 0;

    reconPatchSizes = 0;

    // HR 21.12.18 : A bit of a hack for the moment in order to support the new
    // method (using global point addressing) and old method (syncing across
    // processor BCs) of constructing the shared points. The old method is still
    // needed since global point addressing does not exist when the sharedPoint
    // object is constructed from reconstructParMesh.
    //
    // TODO: Unify the methods by constructing global point addressing from the
    // mesh pieces when. Or even better, calculate procPointAddressing directly
    // which will simplify the mesh merging immensely.
    autoPtr<sharedPoints> sharedDataPtr;
    if(globalPointIndex_.size())
    {
        // Determine globally shared points using global point
        // adressing
        sharedDataPtr.set(new sharedPoints(meshes_, globalPointIndex_));
    }
    else
    {
        // Prepare handling for globally shared points.  This is equivalent
        // to parallel processor points, but working on a PtrList of meshes
        // on the same processor
        sharedDataPtr.set(new sharedPoints(meshes_));
    }
    sharedPoints& sharedData = sharedDataPtr();

    // Record global point index for shared points
    labelList globalPointMapping(sharedData.nGlobalPoints(), -1);

    // Before assembling the meshes, create unique ordering for all passive
    // processor patches that will be merged.  HJ, 6/May/2018

    // This list gives insert order for every passive processor patch in the
    // boundary mesh of every processor mesh.
    // The index gives the location in the reconstructed patch where the
    // face should be inserted
    // For other patch types, the list is empty
    labelListListList passivePatchInsertOrder(meshes_.size());


    // Memory management
    {
        // This list provides insert offset for the patch to unpack the list
        labelListList passivePatchInsertOffset(meshes_.size());

        // This list records all global indices from passive faces
        // (from multiple proc meshes and patches) that will end up together
        // It refers to the new patches
        labelListList globalIndexPerNewProcPatch(reconPatchTypes.size());

        forAll (meshes_, procI)
        {
            if (meshes_.set(procI))
            {
                const polyBoundaryMesh& procPatches =
                    meshes_[procI].boundaryMesh();

                passivePatchInsertOrder[procI].setSize(procPatches.size());
                passivePatchInsertOffset[procI].setSize(procPatches.size());

                forAll (procPatches, patchI)
                {
                    if (isA<passiveProcessorPolyPatch>(procPatches[patchI]))
                    {
                        const passiveProcessorPolyPatch& curPatch =
                            refCast<const passiveProcessorPolyPatch>
                            (
                                procPatches[patchI]
                            );

                        // Find new patch index in reconstructed boundary
                        const label pnIndex = patchNameLookup.find
                            (curPatch.name())();

                        // Record offset
                        passivePatchInsertOffset[procI][patchI] =
                            globalIndexPerNewProcPatch[pnIndex].size();

                        // Set or append the global index
                        globalIndexPerNewProcPatch[pnIndex].append
                        (
                            curPatch.globalFaceIndex()
                        );
                    }
                }
            }
        }

        // Now all indices and offsets are collected.
        // Sort each individual list and get sorted index order
        // Note: this refers to the new patch
        forAll (globalIndexPerNewProcPatch, newPatchI)
        {
            if (!globalIndexPerNewProcPatch[newPatchI].empty())
            {
                SortableList<label> sLabels
                (
                    globalIndexPerNewProcPatch[newPatchI]
                );

                // Note: on sort, indices are done the wrong way around:
                // old index for new location: Re-pack them
                const labelList& indices = sLabels.indices();

                // Re-use the storage for the insertion index
                labelList& invIndices = globalIndexPerNewProcPatch[newPatchI];

                forAll (invIndices, i)
                {
                    invIndices[indices[i]] = i;
                }
            }
        }

        // Unpack the list for insertion
        forAll (meshes_, procI)
        {
            if (meshes_.set(procI))
            {
                const polyBoundaryMesh& procPatches =
                    meshes_[procI].boundaryMesh();

                forAll (procPatches, patchI)
                {
                    if (isA<passiveProcessorPolyPatch>(procPatches[patchI]))
                    {
                        const passiveProcessorPolyPatch& curPatch =
                            refCast<const passiveProcessorPolyPatch>
                            (
                                procPatches[patchI]
                            );

                        // Find new patch index in reconstructed boundary
                        const label pnIndex = patchNameLookup.find
                            (curPatch.name())();

                        // Collect offset
                        passivePatchInsertOrder[procI][patchI] =
                            labelList::subList
                            (
                                globalIndexPerNewProcPatch[pnIndex],
                                curPatch.size(),
                                passivePatchInsertOffset[procI][patchI]
                            );
                    }
                }
            }
        }
    }

    // Dump first valid mesh without checking
    {
        const label fvmId = firstValidMesh();

        if (debug)
        {
            Pout<< "Dump mesh " << fvmId << endl;
        }

        cellOffset[fvmId] = 0;

        const polyMesh& curMesh = meshes_[fvmId];
        labelList& ppAddr = pointProcAddressing_[fvmId];
        labelList& fpAddr = faceProcAddressing_[fvmId];
        labelList& cpAddr = cellProcAddressing_[fvmId];
        labelList& bpAddr = boundaryProcAddressing_[fvmId];

        // Dump all points into the global point list
        // Reconstruct only live points.  HJ, 7/Mar/2011
        const pointField& curPoints = curMesh.points();
        ppAddr.setSize(curPoints.size());

        forAll (curPoints, pointI)
        {
            reconPoints[nReconPoints] = curPoints[pointI];
            ppAddr[pointI] = nReconPoints;
            nReconPoints++;
        }

        // Collect globally shared point labels
        const labelList& curSpAddr = sharedData.sharedPointAddr()[fvmId];
        const labelList& curSpl = sharedData.sharedPointLabels()[fvmId];

        forAll (curSpAddr, splI)
        {
            // From processor 0, mark points without checking
            globalPointMapping[curSpAddr[splI]] = ppAddr[curSpl[splI]];
        }

        // Dump all internal faces into the list
        const faceList& curFaces = curMesh.allFaces();
        const labelList& curOwner = curMesh.faceOwner();
        const labelList& curNeighbour = curMesh.faceNeighbour();
        fpAddr.setSize(curFaces.size());

        for (label faceI = 0; faceI < curMesh.nInternalFaces(); faceI++)
        {
            // Renumber face in new vertices
            face newFace = curFaces[faceI];
            inplaceRenumber(ppAddr, newFace);

            reconFaces[nReconFaces] = newFace;
            reconOwner[nReconFaces] = curOwner[faceI];//+cellOffset[firstValidMesh()];
            reconNeighbour[nReconFaces] = curNeighbour[faceI];//+cellOffset[firstValidMesh()];

            // Face-processor addressing uses offset of 1 and a turning index
            // If the label is negative, it means the global face points
            // in direction opposite to decomposed face.  HJ, 16/Feb/2011
            fpAddr[faceI] = nReconFaces + 1;
            nReconFaces++;
        }

        // Go through all patches.
        // For processor make internal faces
        // For regular patches dump the faces into patch lists
        // Note: processor patches need to be visited in the increasing
        // neighbour processor index.  HJ, 24/May/2018

        const polyBoundaryMesh& procPatches = curMesh.boundaryMesh();

        // Sort processor patches for insertion
        labelList nrbProcIndex(procPatches.size(), -1);

        forAll (procPatches, patchI)
        {
            if (isA<processorPolyPatch>(procPatches[patchI]))
            {
                // Processor patch: faces become internal faces
                const processorPolyPatch& curPatch =
                    refCast<const processorPolyPatch>(procPatches[patchI]);

                nrbProcIndex[patchI] = curPatch.neighbProcNo();
            }
            else if (isA<passiveProcessorPolyPatch>(procPatches[patchI]))
            {
                // For passive processor patch, faces need to be inserted in
                // the increasing global face index
                const passiveProcessorPolyPatch& curPatch =
                    refCast<const passiveProcessorPolyPatch>
                    (
                        procPatches[patchI]
                    );

                nrbProcIndex[patchI] = curPatch.neighbProcNo();
            }
        }

        // Get sorting order. Note: make a copy of indices because
        // sortable list will be deleted
        const labelList procVisitOrder =
            SortableList<label>(nrbProcIndex).indices();

        forAll (procVisitOrder, pvoI)
        {
            const label& patchI = procVisitOrder[pvoI];

            if (isA<processorPolyPatch>(procPatches[patchI]))
            {
                // Processor patch: faces become internal faces
                const processorPolyPatch& curPatch =
                    refCast<const processorPolyPatch>(procPatches[patchI]);

                // Record boundary-processor addressing: unmapped patch
                bpAddr[patchI] = -1;

                const label patchStart = curPatch.start();

                for
                (
                    label faceI = patchStart;
                    faceI < patchStart + curPatch.size();
                    faceI++
                )
                {
                    // Renumber face in new vertices
                    face newFace = curFaces[faceI];
                    inplaceRenumber(ppAddr, newFace);

                    reconFaces[nReconFaces] = newFace;

                    fpAddr[faceI] = nReconFaces + 1;
                    reconOwner[nReconFaces] = curOwner[faceI];//+cellOffset[firstValidMesh()];

                    // For partially completed neighbour, set nbr to -2
                    // for easier debugging
                    reconNeighbour[nReconFaces] = -2;

                    nReconFaces++;
                }
            }
            else if (isA<passiveProcessorPolyPatch>(procPatches[patchI]))
            {
                // For passive processor patch, faces need to be inserted in
                // the increasing global face index
                const passiveProcessorPolyPatch& curPatch =
                    refCast<const passiveProcessorPolyPatch>
                    (
                        procPatches[patchI]
                    );

                // Find processor patch index in reconstructed boundary
                const label pnIndex = patchNameLookup.find(curPatch.name())();

                // Record boundary-processor addressing: mapped patch
                bpAddr[patchI] = pnIndex;

                faceList& curRpFaces = reconPatchFaces[pnIndex];
                labelList& curRpfOwner = reconPatchOwner[pnIndex];
                label& nRpf = reconPatchSizes[pnIndex];

                const labelList& indices =
                    passivePatchInsertOrder[fvmId][patchI];

                const label patchStart = curPatch.start();

                forAll (indices, i)
                {
                    // Location in reconstructed patch where the face is
                    // inserted
                    const label& insertSlot = indices[i];

                    // Calculate face index depending on the ordering
                    const label faceI = patchStart + i;

                    face newFace = curFaces[faceI];
                    inplaceRenumber(ppAddr, newFace);

                    // Insert into correct slot
                    curRpFaces[insertSlot] = newFace;
                    curRpfOwner[insertSlot] = curOwner[faceI];//+cellOffset[firstValidMesh()];

                    // Temporarily record position of face in the patch.
                    // Offset for nInternalFaces will be added in the end
                    // when the complete list of faces is assembled
                    // HJ, 16/Feb/2011
                    fpAddr[faceI] = insertSlot + 1;
                    nRpf++;
                }
            }
            else
            {
                // Regular patch: dump faces into patch face list

                const polyPatch& curPatch = procPatches[patchI];

                // Find processor patch index in reconstructed boundary
                const label pnIndex = patchNameLookup.find(curPatch.name())();

                // Record boundary-processor addressing: mapped patch
                bpAddr[patchI] = pnIndex;

                faceList& curRpFaces = reconPatchFaces[pnIndex];
                labelList& curRpfOwner = reconPatchOwner[pnIndex];
                label& nRpf = reconPatchSizes[pnIndex];

                const label patchStart = curPatch.start();

                for
                (
                    label faceI = patchStart;
                    faceI < patchStart + curPatch.size();
                    faceI++
                )
                {
                    // Renumber face in new vertices
                    face newFace = curFaces[faceI];
                    inplaceRenumber(ppAddr, newFace);

                    curRpFaces[nRpf] = newFace;
                    curRpfOwner[nRpf] = curOwner[faceI];//+cellOffset[firstValidMesh()];

                    // Temporarily record position of face in the patch.
                    // Offset for nInternalFaces will be added in the end
                    // when the complete list of faces is assembled
                    // HJ, 16/Feb/2011
                    fpAddr[faceI] = nRpf + 1;
                    nRpf++;
                }
            }
        }

        // Cell-processor addressing
        forAll (cpAddr, cellI)
        {
            cpAddr[cellI] = cellI; // + cellOffset[firstValidMesh()];
        }

        // Set cell offset for the next mesh. Note: if the next mesh is not
        // valid, the cell offset is propagated to the next one in the processor
        // loop below.
        if (cellOffset.size() > firstValidMesh() + 1)
        {
            cellOffset[firstValidMesh() + 1] =
                meshes_[firstValidMesh()].nCells();
        }
    }


    // Dump all other meshes, merging the processor boundaries
    for (label procI = firstValidMesh() + 1; procI < meshes_.size(); procI++)
    {
        if (!meshes_.set(procI))
        {
            // Next mesh is not valid: simply propagate cell offset
            if (cellOffset.size() > procI + 1)
            {
                cellOffset[procI + 1] = cellOffset[procI];
            }
        }
        else
        {
            // Valid mesh, combine it
            if (debug)
            {
                Pout<< "Dump mesh " << procI
                    << " cell offset: " << cellOffset[procI]
                    << endl;
            }

            const polyMesh& curMesh = meshes_[procI];
            const polyBoundaryMesh& procPatches = curMesh.boundaryMesh();

            labelList& ppAddr = pointProcAddressing_[procI];
            labelList& fpAddr = faceProcAddressing_[procI];
            labelList& cpAddr = cellProcAddressing_[procI];
            labelList& bpAddr = boundaryProcAddressing_[procI];

            // Point mapping

            // Reconstruct only live points.  HJ, 7/Mar/2011
            const pointField& curPoints = curMesh.points();

            // Set ppAddr to -1, to use as point usage indicators
            ppAddr.setSize(curPoints.size());
            ppAddr = -1;

            // Collect globally shared point labels
            const labelList& curSpAddr = sharedData.sharedPointAddr()[procI];
            const labelList& curSpl = sharedData.sharedPointLabels()[procI];

            // Go through globally shared points.  If already marked, mark
            // their point index.  If not, leave for later

            forAll (curSpAddr, splI)
            {
                // From other processors, check if point is already marked
                // If yes, insert the point label into ppAddr
                if (globalPointMapping[curSpAddr[splI]] > -1)
                {
                    // Point is marked.  Record it locally
                    ppAddr[curSpl[splI]] = globalPointMapping[curSpAddr[splI]];
                }
            }

            // Find points already added via processor patches and mark them
            // in ppAddr

            // Collect point-processor addressing for points on
            // processor patches

            // Go through all processor patches.  For neighbour patches, access
            // owner addressing and dump into ppAddr
            forAll (procPatches, patchI)
            {
                if (isA<processorPolyPatch>(procPatches[patchI]))
                {
                    // Processor patch: faces become internal faces
                    const processorPolyPatch& procPatch =
                        refCast<const processorPolyPatch>(procPatches[patchI]);

                    // If patch is a neighbour, its master has already inserted
                    // the points
                    if (procPatch.slave())
                    {
                        const label masterProcID = procPatch.neighbProcNo();

                        // Get the neighbour side patch
                        const processorPolyPatch& masterProcPatch =
                            neighbourProcPatch(procPatch);

                        // Find the addressing of the master side
                        const labelList& masterPpAddr =
                            pointProcAddressing_[masterProcID];

                        // Assemble neighbour mesh point addressing in matching
                        // order by reversing processor patch faces
                        faceList reversedFaces(procPatch.size());

                        forAll (reversedFaces, faceI)
                        {
                            reversedFaces[faceI] =
                                procPatch[faceI].reverseFace();
                        }

                        primitiveFacePatch reversedPatch
                        (
                            reversedFaces,
                            procPatch.points()
                        );

                        // Insert addressing from master side into
                        // local point addressing.  Each face of reversed patch
                        // now matches the master face.
                        // Note: this is done by visiting faces, since
                        // meshPoints are ordered in increasing order.
                        // HJ, 10/Mar/2011

                        forAll (reversedFaces, faceI)
                        {
                            // Current reverse face
                            const face& curRF = reversedFaces[faceI];

                            // Current master face
                            const face& curMF = masterProcPatch[faceI];

                            forAll (curRF, pointI)
                            {
                                if (ppAddr[curRF[pointI]] == -1)
                                {
                                    // Mapping is established
                                    ppAddr[curRF[pointI]] =
                                        masterPpAddr[curMF[pointI]];
                                }
                                else
                                {
                                    // Mapping already exists.  Check it
                                    if
                                    (
                                        ppAddr[curRF[pointI]]
                                     != masterPpAddr[curMF[pointI]]
                                    )
                                    {
                                        FatalErrorIn
                                        (
                                            "autoPtr<fvMesh> "
                                            "processorMeshesReconstructor::"
                                            "reconstructMesh(const Time& db)"
                                        )   << "Loss of proc sync: proc pair: ("
                                            << procPatch.myProcNo()
                                            << " " << procPatch.neighbProcNo()
                                            << ") point addr: "
                                            << ppAddr[curRF[pointI]] << " and "
                                            << masterPpAddr[curMF[pointI]]
                                            << " for point "
                                            << curRF[pointI]
                                            << " at "
                                            << curMesh.points()[curRF[pointI]]
                                            << nl
                                            << abort(FatalError);
                                    }
                                }
                            }
                        }
                    } // End of "is neighbour"
                } // End of "is processor"
            } // End for all patches

            // Dump unmarked points into the global point list
            label nMergedPoints = 0;

            forAll (curPoints, pointI)
            {
                if (ppAddr[pointI] == -1)
                {
                    // Unmerged point
                    reconPoints[nReconPoints] = curPoints[pointI];
                    ppAddr[pointI] = nReconPoints;
                    nReconPoints++;
                }
                else
                {
                    nMergedPoints++;
                }
            }

            // Dump all internal faces into the list
            const faceList& curFaces = curMesh.allFaces();
            const labelList& curOwner = curMesh.faceOwner();
            const labelList& curNeighbour = curMesh.faceNeighbour();
            fpAddr.setSize(curFaces.size());

            // Collect newly marked globally shared point labels
            forAll (curSpAddr, splI)
            {
                if (globalPointMapping[curSpAddr[splI]] < 0)
                {
                    globalPointMapping[curSpAddr[splI]] = ppAddr[curSpl[splI]];
                }
                else
                {
                    // If point is already set, compare sync
                    // to detect merge error.  Debug
                    if
                    (
                        globalPointMapping[curSpAddr[splI]]
                     != ppAddr[curSpl[splI]]
                    )
                    {
                        FatalErrorIn
                        (
                            "autoPtr<fvMesh> "
                            "processorMeshesReconstructor::"
                            "reconstructMesh(const Time& db)"
                        )   << "Loss of global sync: "
                            << globalPointMapping[curSpAddr[splI]] << " and "
                            << ppAddr[curSpl[splI]] << nl
                            << abort(FatalError);
                    }
                }
            }

            for (label faceI = 0; faceI < curMesh.nInternalFaces(); faceI++)
            {
                // Renumber face in new vertices
                face newFace = curFaces[faceI];
                inplaceRenumber(ppAddr, newFace);

                reconFaces[nReconFaces] = newFace;
                reconOwner[nReconFaces] = curOwner[faceI] + cellOffset[procI];
                reconNeighbour[nReconFaces] = curNeighbour[faceI]
                    + cellOffset[procI];
                fpAddr[faceI] = nReconFaces + 1;
                nReconFaces++;
            }

            // Go through all patches.
            // For processor make internal faces
            // Note: processor patches need to be visited in the increasing
            // neighbour processor index.  HJ, 24/May/2018

            // Sort processor patches for insertion
            labelList nrbProcIndex(procPatches.size(), -1);

            forAll (procPatches, patchI)
            {
                if (isA<processorPolyPatch>(procPatches[patchI]))
                {
                    // Processor patch: faces become internal faces
                    const processorPolyPatch& curPatch =
                        refCast<const processorPolyPatch>(procPatches[patchI]);

                    nrbProcIndex[patchI] = curPatch.neighbProcNo();
                }
                else if (isA<passiveProcessorPolyPatch>(procPatches[patchI]))
                {
                    // For passive processor patch, faces need to be inserted in
                    // the increasing global face index
                    const passiveProcessorPolyPatch& curPatch =
                        refCast<const passiveProcessorPolyPatch>
                        (
                            procPatches[patchI]
                        );

                    nrbProcIndex[patchI] = curPatch.neighbProcNo();
                }
            }

            // Get sorting order.  Note make a copy of indices because
            // sortable list will be deleted
            const labelList procVisitOrder =
                SortableList<label>(nrbProcIndex).indices();

            forAll (procVisitOrder, pvoI)
            {
                const label patchI = procVisitOrder[pvoI];

                if (isA<processorPolyPatch>(procPatches[patchI]))
                {
                    // Processor patch: faces become internal faces
                    const processorPolyPatch& curPatch =
                        refCast<const processorPolyPatch>(procPatches[patchI]);

                    // Record boundary-processor addressing: unmapped patch
                    bpAddr[patchI] = -1;

                    const label patchStart = curPatch.start();

                    // If patch is a master, drop the faces and fill the
                    // owner side addressing
                    if (curPatch.master())
                    {
                        for
                        (
                            label faceI = patchStart;
                            faceI < patchStart + curPatch.size();
                            faceI++
                        )
                        {
                            // Renumber face in new vertices
                            face newFace = curFaces[faceI];
                            inplaceRenumber(ppAddr, newFace);

                            reconFaces[nReconFaces] = newFace;
                            fpAddr[faceI] = nReconFaces + 1;
                            reconOwner[nReconFaces] = curOwner[faceI]
                                + cellOffset[procI];

                            // For partially completed neighbour, set nbr to -2
                            // for easier debugging
                            reconNeighbour[nReconFaces] = -2;

                            nReconFaces++;
                        }
                    }
                    else
                    {
                        // Fill the addressing for the neighbour side

                        const label masterProcID = curPatch.neighbProcNo();

                        // Get local face-cell addressing: it will become
                        // a neighbour
                        // addressing for the already inserted faces
                        const labelList procFaceCells = curPatch.faceCells();

                        // Get the neighbour side patch
                        const processorPolyPatch& masterProcPatch =
                            neighbourProcPatch(curPatch);

                        // Find the addressing of the master side
                        // and insert the neighbour with offset

                        const labelList& masterFaceAddr =
                            faceProcAddressing_[masterProcID];

                        for
                        (
                            label faceI = patchStart;
                            faceI < patchStart + curPatch.size();
                            faceI++
                        )
                        {
                            const label faceInPatch = faceI - patchStart;

                            // Calculate master index
                            const label masterIndex =
                                masterProcPatch.start() + faceInPatch;

                            const label masterFp =
                                masterFaceAddr[masterIndex] - 1;

                            // Record face-cells for the neighbour
                            fpAddr[faceI] = -masterFaceAddr[masterIndex];

                            reconNeighbour[masterFp] =
                                procFaceCells[faceInPatch] + cellOffset[procI];
                        }
                    }
                }
                else if (isA<passiveProcessorPolyPatch>(procPatches[patchI]))
                {
                    // For passive processor patch, faces need to be inserted in
                    // the increasing global face index
                    const passiveProcessorPolyPatch& curPatch =
                        refCast<const passiveProcessorPolyPatch>
                        (
                            procPatches[patchI]
                        );

                    // Find processor patch index in reconstructed boundary
                    const label pnIndex =
                        patchNameLookup.find(curPatch.name())();

                    // Record boundary-processor addressing: mapped patch
                    bpAddr[patchI] = pnIndex;

                    faceList& curRpFaces = reconPatchFaces[pnIndex];
                    labelList& curRpfOwner = reconPatchOwner[pnIndex];
                    label& nRpf = reconPatchSizes[pnIndex];

                    const labelList& indices =
                        passivePatchInsertOrder[procI][patchI];

                    const label patchStart = curPatch.start();

                    forAll (indices, i)
                    {
                        // Location in reconstructed patch where the face
                        // is inserted
                        const label insertSlot = indices[i];

                        // Calculate face index depending on the ordering
                        const label faceI = patchStart + i;

                        face newFace = curFaces[faceI];

                        inplaceRenumber(ppAddr, newFace);

                        // Insert into correct slot
                        curRpFaces[insertSlot] = newFace;
                        curRpfOwner[insertSlot] =
                            curOwner[faceI] + cellOffset[procI];

                        // Temporarily record position of face in the patch.
                        // Offset for nInternalFaces will be added in the end
                        // when the complete list of faces is assembled
                        // HJ, 16/Feb/2011
                        fpAddr[faceI] = insertSlot + 1;
                        nRpf++;
                    }
                }
                else
                {
                    // Regular patch: dump faces into patch face list

                    const polyPatch& curPatch = procPatches[patchI];

                    // Find processor patch index in reconstructed boundary
                    const label pnIndex =
                        patchNameLookup.find(curPatch.name())();

                    // Record boundary-processor addressing: mapped patch
                    bpAddr[patchI] = pnIndex;

                    faceList& curRpFaces = reconPatchFaces[pnIndex];
                    labelList& curRpfOwner = reconPatchOwner[pnIndex];
                    label& nRpf = reconPatchSizes[pnIndex];

                    const label patchStart = curPatch.start();

                    for
                    (
                        label faceI = patchStart;
                        faceI < patchStart + curPatch.size();
                        faceI++
                    )
                    {
                        // Renumber face in new vertices
                        face newFace = curFaces[faceI];
                        inplaceRenumber(ppAddr, newFace);

                        curRpFaces[nRpf] = newFace;
                        curRpfOwner[nRpf] = curOwner[faceI] + cellOffset[procI];

                        // Temporarily record position of face in the patch.
                        // Offset for nInternalFaces will be added in the end
                        // when the complete list of faces is assembled
                        // HJ, 16/Feb/2011
                        fpAddr[faceI] = nRpf + 1;
                        nRpf++;
                    }
                }
            }

            // Cell-processor addressing
            forAll (cpAddr, cellI)
            {
                cpAddr[cellI] = cellI + cellOffset[procI];
            }

            // Set cell offset for the next mesh
            if (cellOffset.size() > procI + 1)
            {
                cellOffset[procI + 1] = cellOffset[procI] + meshes_[procI].nCells();
            }
        }
    }

    // Resize the lists
    reconPoints.setSize(nReconPoints);

    // Resize the neighbour list to the size of internalFaces
    label nInternalFaces = nReconFaces;

    reconNeighbour.setSize(nInternalFaces);

    // Resize the patch lists
    forAll (reconPatchFaces, patchI)
    {
        reconPatchFaces[patchI].setSize(reconPatchSizes[patchI]);
        reconPatchOwner[patchI].setSize(reconPatchSizes[patchI]);
    }

    // Complete the global list of faces
    labelList reconPatchStarts(reconPatchSizes, 0);

    // Copy the faces into face list
    forAll (reconPatchFaces, patchI)
    {
        reconPatchStarts[patchI] = nReconFaces;

        const faceList& curPatchFaces = reconPatchFaces[patchI];
        const labelList& curPatchOwner = reconPatchOwner[patchI];

        forAll (curPatchFaces, fI)
        {
            reconFaces[nReconFaces] = curPatchFaces[fI];
            reconOwner[nReconFaces] = curPatchOwner[fI];
            nReconFaces++;
        }
    }

    reconFaces.setSize(nReconFaces);
    reconOwner.setSize(nReconFaces);

    // Mesh assembly completed

    if (!Pstream::parRun())
    {
        Info<< "Global mesh size (final): " << nl
            << "    nPoints = " << reconPoints.size() << nl
            << "    nFaces = " << reconFaces.size() << nl
            << "    nCells = " << nReconCells << nl
            << "    nPatches = " << reconPatchSizes.size() << nl
            << "    nPatchFaces = " << reconPatchSizes << endl;
    }
    else if (debug)
    {
        Pout<< "Local mesh size (final): " << nl
            << "    nPoints = " << reconPoints.size() << nl
            << "    nFaces = " << reconFaces.size() << nl
            << "    nCells = " << nReconCells << nl
            << "    nPatches = " << reconPatchSizes.size() << nl
            << "    nPatchFaces = " << reconPatchSizes << endl;
    }

    // Renumber the face-processor addressing list for all pieces
    // now that the number of internal faces is known
    forAll (meshes_, procI)
    {
        if (meshes_.set(procI))
        {
            // Get processor mesh and boundary
            const polyMesh& curMesh = meshes_[procI];
            const polyBoundaryMesh& procPatches = curMesh.boundaryMesh();

            // Get face-processor addressing for corrent prorcessor
            labelList& fpAddr = faceProcAddressing_[procI];

            const labelList& bpAddr = boundaryProcAddressing_[procI];

            forAll (procPatches, patchI)
            {
                if (!isA<processorPolyPatch>(procPatches[patchI]))
                {
                    // Get master processor patch
                    const label reconPatchID = bpAddr[patchI];

                    // Skip processor patches: bpAddr = -1
                    if (reconPatchID > -1)
                    {
                        const label reconStart = reconPatchStarts[reconPatchID];

                        const polyPatch& curPatch = procPatches[patchI];

                        const label patchStart = curPatch.start();

                        for
                        (
                            label faceI = patchStart;
                            faceI < patchStart + curPatch.size();
                            faceI++
                        )
                        {
                            // Add patch start
                            fpAddr[faceI] += reconStart;
                        }
                    }
                }
            }
        }
    }

    // Create global mesh with given region name
    autoPtr<fvMesh> globalMeshPtr
    (
        new fvMesh
        (
            IOobject
            (
                meshName_,
                db.timeName(),
                db,
                IOobject::NO_READ
            ),
            xferCopy(reconPoints),
            xferCopy(reconFaces),
            xferCopy(reconOwner),
            xferCopy(reconNeighbour),
            true
            // false                       // Do not sync par
        )
    );
    fvMesh& globalMesh = globalMeshPtr();

    // Create patch list by cloning meshes.  If all processors hold all live
    // patches, it is sufficient to rebuilt the patches only from the first
    // valid processor
    // Note:
    List<polyPatch*> reconPatches(nReconPatches, nullptr);

    forAll (meshes_, procI)
    {
        if (meshes_.set(procI))
        {
            const polyBoundaryMesh& procPatches = meshes_[procI].boundaryMesh();

            forAll (procPatches, patchI)
            {
                // Processor patches have disappeared: skip them
                if (!isA<processorPolyPatch>(procPatches[patchI]))
                {
                    // Find processor patch index in reconstructed boundary
                    const label pnIndex = patchNameLookup.find
                    (
                        procPatches[patchI].name()
                    )();

                    // Check if the patch has already been set
                    if (reconPatches[pnIndex] == nullptr)
                    {
                        // Patch not set: clone it
                        // Note: watch indices: setting pnIndex from patchI

                        if (isA<passiveProcessorPolyPatch>(procPatches[patchI]))
                        {
                            // For a passive processor patch, create new
                            // processor patch
                            const passiveProcessorPolyPatch& ppPatch =
                                refCast<const passiveProcessorPolyPatch>
                                (
                                    procPatches[patchI]
                                );

                            reconPatches[pnIndex] = new processorPolyPatch
                            (
                                ppPatch.name(),
                                reconPatchSizes[pnIndex],
                                reconPatchStarts[pnIndex],
                                pnIndex,
                                globalMesh.boundaryMesh(),
                                Pstream::myProcNo(),   // Use correct local proc
                                ppPatch.neighbProcNo()
                            );
                        }
                        else
                        {
                            // Regular patch: clone
                            reconPatches[pnIndex] =
                                procPatches[patchI].clone
                                (
                                    globalMesh.boundaryMesh(),
                                    patchI,
                                    reconPatchSizes[pnIndex],
                                    reconPatchStarts[pnIndex]
                                ).ptr();
                        }
                    }
                }
            }
        }
    }

    // Check the list and fill in the missing slots
    forAll (reconPatches, patchI)
    {
        if (reconPatches[patchI] == nullptr)
        {
            // Patch not set.  Check its type
            FatalErrorIn
            (
                "autoPtr<fvMesh> processorMeshesReconstructor::"
                "reconstructMesh(const Time& db)"
            )   << "Reconstructed patch " << patchI
                << " name " << reconPatchNames[patchI]
                << " type " << reconPatchTypes[patchI]
                << " not set."
                << abort(FatalError);
        }
    }

    // Add boundary patches to polyMesh and fvMesh
    // Note: Mark boundary as invalid to disable analysis
    // due to the presence of old/new patches
    globalMesh.addFvPatches(reconPatches, false);

    // Recombine cell, face and point zones.
    // Note 1: all zones have to be present on same processors in the same
    // order. This is the result of the decomposition. See
    // domainDecomposition::processorMesh member function
    // Note 2: the code below could be written in a generic way by using a
    // template helper member function, but it's not straightforward since we
    // don't have a list of ZoneMeshes for all processors

    // Get index for the first valid mesh
    const label fvmID = firstValidMesh();

    // First pass: count maximum number of cells, faces and points in zones
    labelList nCellsPerZone(meshes_[fvmID].cellZones().size(), 0);
    labelList nFacesPerZone(meshes_[fvmID].faceZones().size(), 0);
    labelList nPointsPerZone(meshes_[fvmID].pointZones().size(), 0);

    forAll (meshes_, procI)
    {
        if (meshes_.set(procI))
        {
            // Grab the current mesh
            const polyMesh& curMesh = meshes_[procI];

            // PART 1: Cell zones. Scope for clarity and safety
            {
                const cellZoneMesh& cz = curMesh.cellZones();

                forAll (cz, zoneI)
                {
                    // Count number of cells in the zone
                    nCellsPerZone[zoneI] += cz[zoneI].size();
                }
            }

            // PART 2: Face zones. Scope for clarity and safety
            {
                const faceZoneMesh& fz = curMesh.faceZones();

                forAll (fz, zoneI)
                {
                    // Count number of faces in the zone
                    nFacesPerZone[zoneI] += fz[zoneI].size();
                }
            }

            // PART 3: Point zones. Scope for clarity and safety
            {
                const pointZoneMesh& pz = curMesh.pointZones();

                forAll (pz, zoneI)
                {
                    // Count number of points in the zone
                    nPointsPerZone[zoneI] += pz[zoneI].size();
                }
            }
        } // End if processor mesh set
    } // End for all processor meshes

    // Second pass: redistribute cells, faces and points in zones

    // Create lists that contain labels of all cells/faces/points in a given
    // zone coming from different processor meshes. Index is the zoneID, which
    // is the same for all processor bits, while the other list collects all the
    // cells/faces/points in a given zone.
    labelListList reconCellZones(nCellsPerZone.size());
    labelListList reconFaceZones(nFacesPerZone.size());
    List<boolList> reconFaceZoneFlips(nFacesPerZone.size());
    labelListList reconPointZones(nPointsPerZone.size());

    // Size the lists appropriatelly for each zone
    forAll (reconCellZones, zoneI)
    {
        reconCellZones[zoneI].setSize(nCellsPerZone[zoneI], -1);
    }
    forAll (reconFaceZones, zoneI)
    {
        reconFaceZones[zoneI].setSize(nFacesPerZone[zoneI], -1);
        reconFaceZoneFlips[zoneI].setSize(nFacesPerZone[zoneI], false);
    }
    forAll (reconPointZones, zoneI)
    {
        reconPointZones[zoneI].setSize(nPointsPerZone[zoneI], -1);
    }

    // Reset counting lists for indexing during list item assignement and for
    // collecting the final number of faces/points in face/pointZones (since
    // these can be actually fewer if e.g. a processor face becomes an internal
    // face).
    nCellsPerZone = 0;
    nFacesPerZone = 0;
    nPointsPerZone = 0;

    // Loop through all the meshes and collect cells/faces/points in the zones
    forAll (meshes_, procI)
    {
        if (meshes_.set(procI))
        {
            // Grab the current mesh
            const polyMesh& curMesh = meshes_[procI];

            // PART 1: Cell zones. Scope for clarity and safety
            {
                const cellZoneMesh& cz = curMesh.cellZones();

                // Get old-to-new cell addressing for this mesh
                const labelList& curCellProcAddr = cellProcAddressing_[procI];

                forAll (cz, zoneI)
                {
                    // Get "new" zone cell index
                    label& nCells = nCellsPerZone[zoneI];

                    // Reference to the new recon zone
                    labelList& zoneReconCells = reconCellZones[zoneI];

                    // Get all the cells in this zone
                    const labelList& zoneCells = cz[zoneI];

                    // Loop through the cells
                    forAll (zoneCells, i)
                    {
                        // Get cell index in the processor mesh
                        const label& oldCellI = zoneCells[i];

                        // Get cell index in the new mesh
                        const label& newCellI = curCellProcAddr[oldCellI];

                        // Redundancy check: check if the the newCellI is
                        // -1. This should not happen because cells have perfect
                        // 1-to-1 mapping
                        if (newCellI != -1)
                        {
                            // Insert the cell in the new recon zone and
                            // increment the counter
                            zoneReconCells[nCells++] = newCellI;
                        }
                        else
                        {
                            FatalErrorIn
                            (
                                "autoPtr<fvMesh>"
                                "\n processorMeshesReconstructor::"
                                "reconstructMesh(const Time& db)"
                            )   << "Found unmapped cell while reconstructing"
                                << " cell zones."
                                << nl
                                << "Cell from processor: " << procI << nl
                                << "Cell zone name: " << cz[zoneI].name() << nl
                                << "Index in the cell zone: " << i << nl
                                << "Old cell index: " << oldCellI << nl
                                << "New cell index: " << newCellI
                                << abort(FatalError);
                        }

                    } // End for all cells in the zone
                } // End for all cell zones
            } // End scope for cell zone handling

            // PART 2: Face zones. Scope for clarity and safety
            {
                const faceZoneMesh& fz = curMesh.faceZones();

                // Get old-to-new face addressing for this mesh
                const labelList& curFaceProcAddr = faceProcAddressing_[procI];

                forAll (fz, zoneI)
                {
                    // Get "new" zone face index
                    label& nFaces = nFacesPerZone[zoneI];

                    // Reference to the new recon zone and flips in the zone
                    labelList& zoneReconFaces = reconFaceZones[zoneI];
                    boolList& zoneReconFaceFlips = reconFaceZoneFlips[zoneI];

                    // Get all the faces in this zone and also their flips
                    const labelList& zoneFaces = fz[zoneI];
                    const boolList& zoneFlipMap = fz[zoneI].flipMap();

                    // Loop through the faces
                    forAll (zoneFaces, i)
                    {
                        // Get the face index in the processor mesh
                        const label& oldFaceI = zoneFaces[i];

                        // Get the face index in the new, reconstructed mesh
                        const label& newFaceI = curFaceProcAddr[oldFaceI];

                        // Check if the face is mapped.
                        // Note:
                        // 1. Need to decrement by 1 because of the face turning
                        // 2. No need to handle negative new indices coming from
                        //    slave processor because we'd end up with
                        //    duplicate entries (two faces on two processors
                        //    merged into a single one)
                        if (newFaceI > 0)
                        {
                            // This is a face that's been correctly
                            // mapped, insert the face in the new zone
                            zoneReconFaces[nFaces] = newFaceI - 1;

                            // Also store the flip map of the face. We don't
                            // need to check whether the flip map has been
                            // preserved because we only get the combined faces
                            // that are inserted from master side.
                            zoneReconFaceFlips[nFaces] = zoneFlipMap[i];

                            // Increment the number of faces for this zone
                            ++nFaces;
                        }
                    } // End for all faces in the zone
                } // End for all face zones
            } // End scope for face zone handling

            // PART 3: Point zones. Scope for clarity and safety
            {
                const pointZoneMesh& fz = curMesh.pointZones();

                // Get old-to-new point addressing for this mesh
                const labelList& curPointProcAddr = pointProcAddressing_[procI];

                forAll (fz, zoneI)
                {
                    // Get "new" zone point index
                    label& nPoints = nPointsPerZone[zoneI];

                    // Reference to the new recon zone
                    labelList& zoneReconPoints = reconPointZones[zoneI];

                    // Get all the points in this zone
                    const labelList& zonePoints = fz[zoneI];

                    // Loop through the points
                    forAll (zonePoints, i)
                    {
                        // Get point index in the processor mesh
                        const label& oldPointI = zonePoints[i];

                        // Get point index in the new mesh
                        const label& newPointI = curPointProcAddr[oldPointI];

                        // Check if the point is mapped
                        if (newPointI != -1)
                        {
                            // Insert the point in the new recon zone and
                            // increment the counter
                            zoneReconPoints[nPoints++] = newPointI;
                        }

                    } // End for all points in the zone
                } // End for all point zones
            } // End scope for point zone handling
        } // End if the processor mesh is set
    } // End for all processor meshes

    // We need to resize the face and point zones to number of inserted
    // faces/points because not all faces and points need to be
    // inserted. There's nothing to do for cell zones because these are always
    // mapped uniquely one-to-one
    forAll (reconFaceZones, zoneI)
    {
        reconFaceZones[zoneI].setSize(nFacesPerZone[zoneI]);
        reconFaceZoneFlips[zoneI].setSize(nFacesPerZone[zoneI]);
    }
    forAll (reconPointZones, zoneI)
    {
        reconPointZones[zoneI].setSize(nPointsPerZone[zoneI]);
    }

    // Now we have all the zones as ordinary lists without possible duplicate
    // faces and points due to merging of processor boundaries. Create zone
    // meshes

    // PART 1: Cell zones
    List<cellZone*> reconCz(reconCellZones.size());

    // Loop through all the cell zones and create them
    forAll (reconCz, zoneI)
    {
        // Notes:
        // 1. Grab the name from the respective zone in the first valid mesh
        // 2. Transfer the list of cell IDs, invalidating reconCellZones[zoneI]
        reconCz[zoneI] = new cellZone
        (
            meshes_[fvmID].cellZones()[zoneI].name(),
            reconCellZones[zoneI].xfer(),
            zoneI,
            globalMesh.cellZones()
        );
    }

    // PART 2: Face zones
    List<faceZone*> reconFz(reconFaceZones.size());

    // Loop through all the face zones and create them
    forAll (reconFz, zoneI)
    {
        // Notes:
        // 1. Grab the name from the respective zone in the first valid mesh
        // 2. Transfer the list of face IDs, invalidating reconFaceZones[zoneI]
        reconFz[zoneI] = new faceZone
        (
            meshes_[fvmID].faceZones()[zoneI].name(),
            reconFaceZones[zoneI].xfer(),
            reconFaceZoneFlips[zoneI].xfer(),
            zoneI,
            globalMesh.faceZones()
        );
    }

    // PART 3: Point zones
    List<pointZone*> reconPz(reconPointZones.size());

    // Loop through all the point zones and create them
    forAll (reconPz, zoneI)
    {
        // Notes:
        // 1. Grab the name from the respective zone in the first valid mesh
        // 2. Transfer the list of point IDs, invalidating reconPointZones[zoneI]
        reconPz[zoneI] = new pointZone
        (
            meshes_[fvmID].pointZones()[zoneI].name(),
            reconPointZones[zoneI].xfer(),
            zoneI,
            globalMesh.pointZones()
        );
    }

    // Add the zones into the mesh
    globalMesh.addZones(reconPz, reconFz, reconCz);

    // All done, return the global mesh pointer
    return globalMeshPtr;
}



// ************************************************************************* //
