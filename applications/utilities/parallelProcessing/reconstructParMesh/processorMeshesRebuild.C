/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "processorMeshesReconstructor.H"
#include "processorPolyPatch.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::processorMeshesReconstructor::readMapping()
{
    // Check for mapping
    Info<< "Check for mesh mapping data for instance "
        << meshes_[0].facesInstance() << ".  ";

    bool readOk = true;

    forAll (meshes_, procI)
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

        Info<< "Addressing from files: " << nl;
        forAll (meshes_, procI)
        {
            Info<< "Proc " << procI
                << " point addr: " << pointProcAddressing_[procI].size()
                << " face addr: " << faceProcAddressing_[procI].size()
                << " cell addr: " << cellProcAddressing_[procI].size()
                << " boundary addr: " << boundaryProcAddressing_[procI].size()
                << endl;
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
        pointProcAddressing_[procI].write();
        faceProcAddressing_[procI].write();
        cellProcAddressing_[procI].write();
        boundaryProcAddressing_[procI].write();
    }
}


const Foam::processorPolyPatch&
Foam::processorMeshesReconstructor::neighbourProcPatch
(
    const processorPolyPatch& procPatch
) const
{
    const label masterProcID = procPatch.neighbProcNo();

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
                // Found matching patch
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
        // Reconstruct only live points.  HJ, 7/Mar/2011
        const vectorField& procPoints = meshes_[procI].points();

        // Set the cell values in the reconstructed field

        const labelList& pointProcAddressingI = pointProcAddressing()[procI];

        if (pointProcAddressingI.size() != procPoints.size())
        {
            FatalErrorIn("processorMeshes")
                << "problem :"
                << " pointProcAddressingI:" << pointProcAddressingI.size()
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

    mesh.movePoints(newPoints);
    mesh.write();
}


Foam::autoPtr<Foam::fvMesh>
Foam::processorMeshesReconstructor::reconstructMesh(const Time& db)
{
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
                meshes_[0].meshSubDir,
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

    // Prepare patch reconstruction from processor 0
    wordList reconPatchNames;
    wordList reconPatchTypes;
    label nReconPatches = 0;

    // Grab names and types of patches, excluding processor patches
    // Order must be identical on all processors
    {
        const polyBoundaryMesh& proc0Boundary = meshes_[0].boundaryMesh();

        reconPatchNames.setSize(proc0Boundary.size());
        reconPatchTypes.setSize(proc0Boundary.size());
        nReconPatches = 0;

        forAll (proc0Boundary, patchI)
        {
            if (!isA<processorPolyPatch>(proc0Boundary[patchI]))
            {
                // Add patch name to list of patches
                reconPatchNames[nReconPatches] = proc0Boundary[patchI].name();
                reconPatchTypes[nReconPatches] = proc0Boundary[patchI].type();
                nReconPatches++;
            }
        }

        reconPatchNames.setSize(nReconPatches);
        reconPatchTypes.setSize(nReconPatches);
    }

    // Prepare point, face and patch reconstruction
    label nReconPoints = 0;
    label nReconFaces = 0;
    label nReconCells = 0;
    labelList reconPatchSizes(reconPatchTypes.size(), 0);

    forAll (meshes_, procI)
    {
        // Count total number of points and faces
        nReconPoints += meshes_[procI].nPoints();
        nReconFaces += meshes_[procI].allFaces().size();
        nReconCells += meshes_[procI].nCells();

        const polyBoundaryMesh& procPatches = meshes_[procI].boundaryMesh();

        nReconPatches = 0;

        forAll (procPatches, patchI)
        {
            if (!isA<processorPolyPatch>(procPatches[patchI]))
            {
                // Check patch name and type
                if
                (
                    procPatches[patchI].name()
                 != reconPatchNames[nReconPatches]
                 || procPatches[patchI].type()
                 != reconPatchTypes[nReconPatches]
                )
                {
                    FatalErrorIn
                    (
                        "autoPtr<fvMesh> processorMeshesReconstructor::"
                        "reconstructMesh(const Time& db)"
                    )   << "Patch names, types or ordering does not match "
                        << "across processors"
                        << abort(FatalError);
                }

                // Record number of faces in patch
                reconPatchSizes[nReconPatches] += procPatches[patchI].size();
                nReconPatches++;
            }
        }
    }

    Info<< "Estimated max global mesh size (with duplicates): " << nl
        << "    nPoints = " << nReconPoints << nl
        << "    nFaces = " << nReconFaces << nl
        << "    nCells = " << nReconCells << nl
        << "    nPatches = " << nReconPatches << nl
        << "    nPatchFaces = " << reconPatchSizes << endl;

    // Note: for easier debugging, set owner and neighbour to -1
    pointField reconPoints(nReconPoints);
    faceList reconFaces(nReconFaces);
    labelList cellOffset(meshes_.size(), 0);
    labelList reconOwner(nReconFaces, -1);
    labelList reconNeighbour(nReconFaces, -1);
    faceListList reconPatchFaces(reconPatchTypes.size());
    labelListList reconPatchOwner(reconPatchTypes.size());

    forAll (reconPatchFaces, patchI)
    {
        reconPatchFaces[patchI].setSize(reconPatchSizes[patchI]);

        reconPatchOwner[patchI].setSize(reconPatchSizes[patchI]);
        reconPatchOwner[patchI] = -1;
    }

    // Size the mapping arrays
    pointProcAddressing_.setSize(meshes_.size());
    faceProcAddressing_.setSize(meshes_.size());
    cellProcAddressing_.setSize(meshes_.size());
    boundaryProcAddressing_.setSize(meshes_.size());

    // Allocate addressing arrays on all meshes
    forAll (meshes_, procI)
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

    // Reset the counters
    nReconPoints = 0;
    nReconFaces = 0;

    reconPatchSizes = 0;

    // Prepare handling for globally shared points
    labelList globalPointMapping;

    // Dump mesh zero without checking
    {
        cellOffset[0] = 0;

        const polyMesh& curMesh = meshes_[0];
        labelList& ppAddr = pointProcAddressing_[0];
        labelList& fpAddr = faceProcAddressing_[0];
        labelList& cpAddr = cellProcAddressing_[0];
        labelList& bpAddr = boundaryProcAddressing_[0];

        // Prepare handling for global mesh data: set to -1
        globalPointMapping.setSize(curMesh.globalData().nGlobalPoints());
        globalPointMapping = -1;

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
        const labelList& curSpl = curMesh.globalData().sharedPointLabels();

        forAll (curSpl, splI)
        {
            // From processor 0, mark points without checking
            globalPointMapping[curSpl[splI]] = ppAddr[curSpl[splI]];
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
            reconOwner[nReconFaces] = curOwner[faceI];//+cellOffset[0];
            reconNeighbour[nReconFaces] = curNeighbour[faceI];//+cellOffset[0];

            // Face-processor addressing uses offset of 1 and a turning index
            // If the label is negative, it means the global face points
            // in direction opposite to decomposed face.  HJ, 16/Feb/2011
            fpAddr[faceI] = nReconFaces + 1;
            nReconFaces++;
        }

        // Go through all patches.  For regular patches
        // dump the faces into patch lists
        const polyBoundaryMesh& procPatches = curMesh.boundaryMesh();

        forAll (procPatches, patchI)
        {
            if (isA<processorPolyPatch>(procPatches[patchI]))
            {
                // Processor patch: faces become internal faces
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(procPatches[patchI]);

                for
                (
                    label faceI = procPatch.start();
                    faceI < procPatch.start() + procPatch.size();
                    faceI++
                )
                {
                    // Renumber face in new vertices
                    face newFace = curFaces[faceI];
                    inplaceRenumber(ppAddr, newFace);

                    reconFaces[nReconFaces] = newFace;

                    fpAddr[faceI] = nReconFaces + 1;
                    reconOwner[nReconFaces] = curOwner[faceI];//+cellOffset[0];

                    // For partially completed neighbour, set nbr to -2
                    // for easier debugging
                    reconNeighbour[nReconFaces] = -2;

                    nReconFaces++;
                }
            }
            else
            {
                // Regular patch: dump faces into patch face list
                faceList& curRpFaces = reconPatchFaces[patchI];
                labelList& curRpfOwner = reconPatchOwner[patchI];
                label& nRpf = reconPatchSizes[patchI];

                const polyPatch& curPatch = procPatches[patchI];

                for
                (
                    label faceI = curPatch.start();
                    faceI < curPatch.start() + curPatch.size();
                    faceI++
                )
                {
                    // Renumber face in new vertices
                    face newFace = curFaces[faceI];
                    inplaceRenumber(ppAddr, newFace);

                    curRpFaces[nRpf] = newFace;
                    curRpfOwner[nRpf] = curOwner[faceI];//+cellOffset[0];

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
            cpAddr[cellI] = cellI; // + cellOffset[0];
        }

        // Sort out boundary addressing: i for live patches, -1 for processor
        bpAddr = -1;

        // Note: loop over mapped patches
        forAll (reconPatchSizes, patchI)
        {
            bpAddr[patchI] = patchI;
        }
    }


    // Dump all other meshes, merging the processor boundaries

    for (label procI = 1; procI < meshes_.size(); procI++)
    {
        // Grab cell offset from previous offset and mesh size
        cellOffset[procI] =
            cellOffset[procI - 1] + meshes_[procI - 1].nCells();

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

        // Find points already added via processor patches and mark them
        // in ppAddr

        // Collect point-processor addressing for points on processor patches

        // Go through all patches.  For neighbour patches, access
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
                if (procPatch.neighbour())
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
                        reversedFaces[faceI] = procPatch[faceI].reverseFace();
                    }

                    primitiveFacePatch reversedPatch
                    (
                        reversedFaces,
                        procPatch.points()
                    );

                    // Insert addressing from master side into
                    // local point addressing.  Each face of reversed patch
                    // now matches the master face.
                    // Note: this is done by visiting faces, since meshPoints
                    // are ordered in increasing order.  HJ, 10/Mar/2011

                    forAll (reversedFaces, faceI)
                    {
                        // Current reverse face
                        const face& curRF = reversedFaces[faceI];

                        // Current master face
                        const face& curMF = masterProcPatch[faceI];

                        forAll (curRF, pointI)
                        {
                            // Mapping is established
                            ppAddr[curRF[pointI]] =
                                masterPpAddr[curMF[pointI]];
                        }
                    }
                } // End of "is neighbour"
            } // End of "is processor"
        }

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

        Info<< "Processor " << procI << " merged " << nMergedPoints
            << " points out of local " << curPoints.size()
            << " and total " << nReconPoints << endl;

        // Dump all internal faces into the list
        const faceList& curFaces = curMesh.allFaces();
        const labelList& curOwner = curMesh.faceOwner();
        const labelList& curNeighbour = curMesh.faceNeighbour();
        fpAddr.setSize(curFaces.size());

        // Collect globally shared point labels
        const labelList& curSpl = curMesh.globalData().sharedPointLabels();

        forAll (curSpl, splI)
        {
            // From other processors, check if point is already marked
            // If not, mark it; otherwise compare (and correct?) with local
            // mark
            if (globalPointMapping[curSpl[splI]] < 0)
            {
                globalPointMapping[curSpl[splI]] = ppAddr[curSpl[splI]];
            }
            else
            {
                // Compare.  Is this needed - should always be OK.
                if (globalPointMapping[curSpl[splI]] != ppAddr[curSpl[splI]])
                {
                    WarningIn
                    (
                        "autoPtr<fvMesh> "
                        "processorMeshesReconstructor::"
                        "reconstructMesh(const Time& db)"
                    )   << "Loss of sync???"
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

        // Go through all patches.  For regular patches
        // dump the faces into patch lists
        forAll (procPatches, patchI)
        {
            if (isA<processorPolyPatch>(procPatches[patchI]))
            {
                // Processor patch: faces become internal faces
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(procPatches[patchI]);

                // If patch is a master, drop the faces and fill the
                // owner side addressing
                if (procPatch.owner())
                {
                    for
                    (
                        label faceI = procPatch.start();
                        faceI < procPatch.start() + procPatch.size();
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

                    const label masterProcID = procPatch.neighbProcNo();

                    // Get local face-cell addressing: it will become neighbour
                    // addressing for the already inserted faces
                    const labelList procFaceCells = procPatch.faceCells();

                    // Get the neighbour side patch
                    const processorPolyPatch& masterProcPatch =
                        neighbourProcPatch(procPatch);

                    // Find the addressing of the master side
                    // and insert the neighbour with offset

                    const labelList& masterFaceAddr =
                        faceProcAddressing_[masterProcID];

                    for
                    (
                        label faceI = procPatch.start();
                        faceI < procPatch.start() + procPatch.size();
                        faceI++
                    )
                    {
                        label faceInPatch = faceI - procPatch.start();

                        // Calculate master index
                        label masterIndex = masterProcPatch.start()
                            + faceInPatch;

                        label masterFp = masterFaceAddr[masterIndex] - 1;

                        // Record face-cells for the neighbour
                        fpAddr[faceI] = -masterFaceAddr[masterIndex];

                        reconNeighbour[masterFp] =
                            procFaceCells[faceInPatch] + cellOffset[procI];
                    }
                }
            }
            else
            {
                // Regular patch: dump faces into patch face list
                faceList& curRpFaces = reconPatchFaces[patchI];
                labelList& curRpfOwner = reconPatchOwner[patchI];
                label& nRpf = reconPatchSizes[patchI];

                const polyPatch& curPatch = procPatches[patchI];

                for
                (
                    label faceI = curPatch.start();
                    faceI < curPatch.start() + curPatch.size();
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

        // Sort out boundary addressing: i for live patches, -1 for processor
        bpAddr = -1;

        // Note: loop over mapped patches
        forAll (reconPatchSizes, patchI)
        {
            bpAddr[patchI] = patchI;
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

    Info<< "Global mesh size (final): " << nl
        << "    nPoints = " << reconPoints.size() << nl
        << "    nFaces = " << reconFaces.size() << nl
        << "    nCells = " << nReconCells << nl
        << "    nPatches = " << reconPatchSizes.size() << nl
        << "    nPatchFaces = " << reconPatchSizes << endl;

    // Renumber the face-processor addressing list for all pieces
    // now that the number of internal faces is known
    forAll (meshes_, procI)
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

                    for
                    (
                        label faceI = curPatch.start();
                        faceI < curPatch.start() + curPatch.size();
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
            xferCopy(reconNeighbour)
        )
    );
    fvMesh& globalMesh = globalMeshPtr();

    // Create patch list using mesh from processor 0
    List<polyPatch*> reconPatches(nReconPatches);

    {
        const polyBoundaryMesh& procPatches = meshes_[0].boundaryMesh();

        reconPatches.setSize(reconPatchSizes.size());

        forAll (reconPatchSizes, patchI)
        {
            reconPatches[patchI] =
                procPatches[patchI].clone
                (
                    globalMesh.boundaryMesh(),
                    patchI,
                    reconPatchSizes[patchI],
                    reconPatchStarts[patchI]
                ).ptr();
        }
    }

    // Add both poly and fv boundary patches
    globalMesh.addFvPatches(reconPatches);

    // TODO: point, face and cell zones

    Info<< "Reconstructed addressing: " << nl;
    forAll (meshes_, procI)
    {
        Info<< "Proc " << procI
            << " point addr: " << pointProcAddressing_[procI].size()
            << " face addr: " << faceProcAddressing_[procI].size()
            << " cell addr: " << cellProcAddressing_[procI].size()
            << " boundary addr: " << boundaryProcAddressing_[procI].size()
            << endl;
    }

    return globalMeshPtr;
}



// ************************************************************************* //
