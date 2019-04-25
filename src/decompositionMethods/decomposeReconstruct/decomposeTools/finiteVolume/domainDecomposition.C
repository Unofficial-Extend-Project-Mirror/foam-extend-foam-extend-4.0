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

#include "domainDecomposition.H"
#include "foamTime.H"
#include "dictionary.H"
#include "labelIOList.H"
#include "labelIOField.H"
#include "processorPolyPatch.H"
#include "passiveProcessorPolyPatch.H"
#include "fvMesh.H"
#include "OSspecific.H"
#include "Map.H"
#include "DynamicList.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::domainDecomposition, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::domainDecomposition::domainDecomposition
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    decompositionDict_(dict),
    nProcs_(readInt(decompositionDict_.lookup("numberOfSubdomains"))),
    distributed_(false),
    gfIndex_(mesh_),
    gpIndex_(mesh_),
    cellToProc_(mesh_.nCells()),
    patchNbrCellToProc_(mesh_.boundaryMesh().size()),
    patchNbrFaceCells_(mesh_.boundaryMesh().size()),
    procPointAddressing_(nProcs_),
    procFaceAddressing_(nProcs_),
    nInternalProcFaces_(nProcs_),
    nLiveProcFaces_(nProcs_),
    procCellAddressing_(nProcs_),
    procBoundaryAddressing_(nProcs_),
    procPatchSize_(nProcs_),
    procPatchStartIndex_(nProcs_),
    procNeighbourProcessors_(nProcs_),
    procProcessorPatchSize_(nProcs_),
    procProcessorPatchStartIndex_(nProcs_),
    globallySharedPoints_(0),
    cyclicParallel_(false)
{
    if (decompositionDict_.found("distributed"))
    {
        distributed_ = Switch(decompositionDict_.lookup("distributed"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::domainDecomposition::~domainDecomposition()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fvMesh> Foam::domainDecomposition::processorMesh
(
    const label procI,
    const Time& processorDb,
    const word& regionName,
    const bool createPassiveProcPatches
) const
{
    // Create processor points
    const labelList& curPointLabels = procPointAddressing_[procI];

    // Access list of all points in the mesh.  HJ, 27/Mar/2009
    const pointField& meshPoints = mesh_.allPoints();

    labelList pointLookup(meshPoints.size(), -1);

    pointField procPoints(curPointLabels.size());

    forAll (curPointLabels, pointi)
    {
        procPoints[pointi] = meshPoints[curPointLabels[pointi]];

        pointLookup[curPointLabels[pointi]] = pointi;
    }

    // Create processor faces
    const labelList& curFaceLabels = procFaceAddressing_[procI];

    // Access list of all faces in the mesh.  HJ, 27/Mar/2009
    const faceList& meshFaces = mesh_.allFaces();

    labelList faceLookup(meshFaces.size(), -1);

    faceList procFaces(curFaceLabels.size());

    forAll (curFaceLabels, faceI)
    {
        // Mark the original face as used
        // Remember to decrement the index by one (turning index) (see
        // procFaceAddressing_ data member comment). HJ, 5/Dec/2001
        label curF = mag(curFaceLabels[faceI]) - 1;

        faceLookup[curF] = faceI;

        // get the original face
        labelList origFaceLabels;

        if (curFaceLabels[faceI] >= 0)
        {
            // face not turned
            origFaceLabels = meshFaces[curF];
        }
        else
        {
            origFaceLabels = meshFaces[curF].reverseFace();
        }

        // translate face labels into local point list
        face& procFaceLabels = procFaces[faceI];

        procFaceLabels.setSize(origFaceLabels.size());

        forAll (origFaceLabels, pointi)
        {
            procFaceLabels[pointi] = pointLookup[origFaceLabels[pointi]];
        }
    }

    // Create cell lookup
    labelList cellLookup(mesh_.nCells(), -1);
    const labelList& curCellLabels = procCellAddressing_[procI];

    forAll (curCellLabels, cellI)
    {
        cellLookup[curCellLabels[cellI]] = cellI;
    }

    // Get complete owner-neighour addressing in the mesh
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    // Calculate owner and neighbour list
    // Owner list is sized to number of live faces
    // Neighbour list is sized to number of internal faces

    labelList procOwner(nLiveProcFaces_[procI]);

    // Note: loop over owner, not all faces: sizes are different
    forAll (procOwner, faceI)
    {
        // Remember to decrement the index by one (turning index) (see
        // procFaceAddressing_ data member comment). HJ, 5/Dec/2001
        label curF = mag(curFaceLabels[faceI]) - 1;

        if (curFaceLabels[faceI] >= 0)
        {
            procOwner[faceI] = cellLookup[own[curF]];
        }
        else
        {
            procOwner[faceI] = cellLookup[nei[curF]];
        }
    }

    labelList procNeighbour(nInternalProcFaces_[procI]);

    // Note: loop over neighbour, not all faces: sizes are different
    forAll (procNeighbour, faceI)
    {
        // Remember to decrement the index by one (turning index) (see
        // procFaceAddressing_ data member comment). HJ, 5/Dec/2001
        label curF = mag(curFaceLabels[faceI]) - 1;

        if (curFaceLabels[faceI] >= 0)
        {
            procNeighbour[faceI] = cellLookup[nei[curF]];
        }
        else
        {
            procNeighbour[faceI] = cellLookup[own[curF]];
        }
    }

    // Create processor mesh without a boundary
    autoPtr<fvMesh> procMeshPtr
    (
        new fvMesh
        (
            IOobject
            (
                regionName,
                mesh_.pointsInstance(),
                processorDb
            ),
            xferMove(procPoints),
            xferMove(procFaces),
            xferMove(procOwner),
            xferMove(procNeighbour),
            false          // Do not sync par
        )
    );
    fvMesh& procMesh = procMeshPtr();

    // Create processor boundary patches
    const labelList& curBoundaryAddressing = procBoundaryAddressing_[procI];

    const labelList& curPatchSizes = procPatchSize_[procI];

    const labelList& curPatchStarts = procPatchStartIndex_[procI];

    const labelList& curNeighbourProcessors =
        procNeighbourProcessors_[procI];

    const labelList& curProcessorPatchSizes =
        procProcessorPatchSize_[procI];

    const labelList& curProcessorPatchStarts =
        procProcessorPatchStartIndex_[procI];

    const polyPatchList& meshPatches = mesh_.boundaryMesh();

    List<polyPatch*> procPatches
    (
        curPatchSizes.size()
      + curProcessorPatchSizes.size(),
        reinterpret_cast<polyPatch*>(0)
    );

    label nPatches = 0;

    forAll (curPatchSizes, patchi)
    {
        procPatches[nPatches] =
            meshPatches[curBoundaryAddressing[patchi]].clone
            (
                procMesh.boundaryMesh(),
                nPatches,
                curPatchSizes[patchi],
                curPatchStarts[patchi]
            ).ptr();

        nPatches++;
    }

    if (createPassiveProcPatches)
    {
        // Creation of passiveProcessor patches requires a global face index

        // Get global index
        const labelList& globalIndex = gfIndex_.globalLabel();

        forAll (curProcessorPatchSizes, procPatchI)
        {
            // Assemble the global face index for all passive processor patch
            // faces
            labelList patchGlobalIndex(curProcessorPatchSizes[procPatchI]);

            const label curPatchStart = curProcessorPatchStarts[procPatchI];

            forAll (patchGlobalIndex, fI)
            {
                // Remember to decrement the index by one (turning index) (see
                // procFaceAddressing_ data member comment). HJ, 5/Dec/2001
                patchGlobalIndex[fI] =
                    globalIndex[mag(curFaceLabels[curPatchStart + fI]) - 1];
            }

            procPatches[nPatches] =
                new passiveProcessorPolyPatch
                (
                    word("passiveProcBoundary") + Foam::name(procI)
                  + word("to")
                  + Foam::name(curNeighbourProcessors[procPatchI]),
                    curProcessorPatchSizes[procPatchI],
                    curProcessorPatchStarts[procPatchI],
                    nPatches,
                    procMesh.boundaryMesh(),
                    procI,
                    curNeighbourProcessors[procPatchI],
                    patchGlobalIndex
                );

            nPatches++;
        }
    }
    else
    {
        forAll (curProcessorPatchSizes, procPatchI)
        {
            procPatches[nPatches] =
                new processorPolyPatch
                (
                    word("procBoundary") + Foam::name(procI)
                  + word("to")
                  + Foam::name(curNeighbourProcessors[procPatchI]),
                    curProcessorPatchSizes[procPatchI],
                    curProcessorPatchStarts[procPatchI],
                    nPatches,
                    procMesh.boundaryMesh(),
                    procI,
                    curNeighbourProcessors[procPatchI]
                );

            nPatches++;
        }
    }


    // Add boundary patches to polyMesh and fvMesh
    // Note: Mark boundary as invalid to disable analysis
    // due to the presence of old/new patches
    procMesh.addFvPatches(procPatches, false);

    // Create and add zones

    // Note:
    // This coding was all wrong, as each point/face/cell may only
    //  belong to a single zone.
    // Additionally, ordering of points/faces/cells in the processor mesh
    // needs to match the ordering in global mesh zones.  Full rewrite.
    // HJ, 30/Mar/2009

    // Create zones if needed
    if
    (
        !mesh_.pointZones().empty()
     || !mesh_.faceZones().empty()
     || !mesh_.cellZones().empty()
    )
    {
        // Point zones
        List<pointZone*> procPz(mesh_.pointZones().size());

        {
            const pointZoneMesh& pz = mesh_.pointZones();

            // Go through all the zoned points and find out if they
            // belong to a processor.  If so, add it to the zone as
            // necessary
            forAll (pz, zoneI)
            {
                const labelList& zonePoints = pz[zoneI];

                labelList procZonePoints(zonePoints.size());
                label nZonePoints = 0;

                forAll (zonePoints, pointI)
                {
                    const label localIndex =
                        pointLookup[zonePoints[pointI]];

                    if (localIndex >= 0)
                    {
                        // Point live on processor: add to zone
                        procZonePoints[nZonePoints] = localIndex;
                        nZonePoints++;
                    }
                }

                // Add the zone
                procZonePoints.setSize(nZonePoints);

                procPz[zoneI] = new pointZone
                (
                    pz[zoneI].name(),
                    procZonePoints,
                    zoneI,
                    procMesh.pointZones()
                );
            }
        }


        // Face zones
        List<faceZone*> procFz(mesh_.faceZones().size());

        {
            const faceZoneMesh& fz = mesh_.faceZones();

            forAll (fz, zoneI)
            {
                const labelList& zoneFaces = fz[zoneI];
                const boolList& flipMap = fz[zoneI].flipMap();

                // Go through all the zoned faces and find out if they
                // belong to a processor.  If so, add it to the zone as
                // necessary

                labelList procZoneFaces(zoneFaces.size());
                boolList procZoneFaceFlips(zoneFaces.size());
                label nZoneFaces = 0;

                forAll (zoneFaces, faceI)
                {
                    const label localIndex = faceLookup[zoneFaces[faceI]];

                    if (localIndex >= 0)
                    {
                        // Face is present on the processor

                        // Add the face to the zone
                        procZoneFaces[nZoneFaces] = localIndex;

                        // Grab the flip
                        bool flip = flipMap[faceI];

                        if (curFaceLabels[localIndex] < 0)
                        {
                            flip = !flip;
                        }

                        procZoneFaceFlips[nZoneFaces] = flip;
                        nZoneFaces++;
                    }
                }

                // Add the zone
                procZoneFaces.setSize(nZoneFaces);
                procZoneFaceFlips.setSize(nZoneFaces);

                procFz[zoneI] = new faceZone
                (
                    fz[zoneI].name(),
                    procZoneFaces,
                    procZoneFaceFlips,
                    zoneI,
                    procMesh.faceZones()
                );
            }
        }

        // Cell zones
        List<cellZone*> procCz(mesh_.cellZones().size());

        {
            const cellZoneMesh& cz = mesh_.cellZones();

            // Go through all the zoned cells and find out if they
            // belong to a processor.  If so, add it to the zone as
            // necessary

            forAll (cz, zoneI)
            {
                const labelList& zoneCells = cz[zoneI];

                labelList procZoneCells(zoneCells.size());
                label nZoneCells = 0;

                forAll (zoneCells, cellI)
                {
                    const label localIndex = cellLookup[zoneCells[cellI]];

                    if (localIndex >= 0)
                    {
                        procZoneCells[nZoneCells] = localIndex;
                        nZoneCells++;
                    }
                }

                // Add the zone
                procZoneCells.setSize(nZoneCells);

                procCz[zoneI] = new cellZone
                (
                    cz[zoneI].name(),
                    procZoneCells,
                    zoneI,
                    procMesh.cellZones()
                );
            }
        }

        // Add zones
        procMesh.addZones(procPz, procFz, procCz);
    }

    // Return mesh
    return procMeshPtr;
}


Foam::labelList Foam::domainDecomposition::globalPointIndex
(
    const label procI
) const
{
    const labelList& gppi = gpIndex_.globalLabel();
    const labelList& ppAddr = procPointAddressing_[procI];
    labelList globalPointIndex(ppAddr.size());

    forAll (globalPointIndex, gpI)
    {
        globalPointIndex[gpI] = gppi[ppAddr[gpI]];
    }

    return globalPointIndex;
}


bool Foam::domainDecomposition::writeDecomposition()
{
    Info<< "\nConstructing processor meshes" << endl;

    // Make a lookup map for globally shared points
    Map<label> sharedPointLookup(2*globallySharedPoints_.size());

    forAll (globallySharedPoints_, pointi)
    {
        sharedPointLookup.insert(globallySharedPoints_[pointi], pointi);
    }


    // Mark point/faces/cells that are in zones.  Bad coding - removed
    // HJ, 31/Mar/2009

    label totProcFaces = 0;
    label maxProcPatches = 0;
    label maxProcFaces = 0;


    // Note: get cellLevel and pointLevel. Avoid checking whether they exist or
    // not by hand. If they don't exist, simply assume that the level is 0
    const labelIOField globalCellLevel
    (
        IOobject
        (
            "cellLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        labelField(mesh_.nCells(), 0)
    );

    const labelIOField globalPointLevel
    (
        IOobject
        (
            "pointLevel",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        labelField(mesh_.nPoints(), 0)
    );


    // Write out the meshes
    for (label procI = 0; procI < nProcs_; procI++)
    {
        fileName processorCasePath
        (
            mesh_.time().caseName()/fileName(word("processor") + name(procI))
        );

        // make the processor directory
        mkDir(mesh_.time().rootPath()/processorCasePath);

        // create a database
        Time processorDb
        (
            Time::controlDictName,
            mesh_.time().rootPath(),
            processorCasePath,
            "system",
            "constant",
            true
        );

        // Set the precision of the points data to 10
        IOstream::defaultPrecision(10);

        autoPtr<fvMesh> procMeshPtr = processorMesh
        (
            procI,
            processorDb,
            mesh_.polyMesh::name()     // region name of undecomposed mesh
        );
        fvMesh& procMesh = procMeshPtr();

        procMesh.write();

        Info<< endl
            << "Processor " << procI << nl
            << "    Number of cells = " << procMesh.nCells()
            << endl;

        label nBoundaryFaces = 0;
        label nProcPatches = 0;
        label nProcFaces = 0;

        forAll (procMesh.boundaryMesh(), patchi)
        {
            if
            (
                procMesh.boundaryMesh()[patchi].type()
             == processorPolyPatch::typeName
            )
            {
                const processorPolyPatch& ppp =
                refCast<const processorPolyPatch>
                (
                    procMesh.boundaryMesh()[patchi]
                );

                Info<< "    Number of faces shared with processor "
                    << ppp.neighbProcNo() << " = " << ppp.size() << endl;

                nProcPatches++;
                nProcFaces += ppp.size();
            }
            else
            {
                nBoundaryFaces += procMesh.boundaryMesh()[patchi].size();
            }
        }

        Info<< "    Number of processor patches = " << nProcPatches << nl
            << "    Number of processor faces = " << nProcFaces << nl
            << "    Number of boundary faces = " << nBoundaryFaces << endl;

        totProcFaces += nProcFaces;
        maxProcPatches = max(maxProcPatches, nProcPatches);
        maxProcFaces = max(maxProcFaces, nProcFaces);

        // create and write the addressing information
        labelIOList pointProcAddressing
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
            procPointAddressing_[procI]
        );
        pointProcAddressing.write();

        labelIOList faceProcAddressing
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
            procFaceAddressing_[procI]
        );
        faceProcAddressing.write();

        labelIOList cellProcAddressing
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
            procCellAddressing_[procI]
        );
        cellProcAddressing.write();

        labelIOList boundaryProcAddressing
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
            procBoundaryAddressing_[procI]
        );
        boundaryProcAddressing.write();

        // Create and write cellLevel and pointLevel information
        const unallocLabelList& cellMap = cellProcAddressing;
        labelIOField procCellLevel
        (
            IOobject
            (
                "cellLevel",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            labelField(globalCellLevel, cellMap)
        );
        procCellLevel.write();

        const unallocLabelList& pointMap = pointProcAddressing;
        labelIOField procPointLevel
        (
            IOobject
            (
                "pointLevel",
                procMesh.facesInstance(),
                procMesh.meshSubDir,
                procMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            labelField(globalPointLevel, pointMap)
        );
        procPointLevel.write();
    }

    Info<< nl
        << "Number of processor faces = " << totProcFaces/2 << nl
        << "Max number of processor patches = " << maxProcPatches << nl
        << "Max number of faces between processors = " << maxProcFaces
        << endl;

    return true;
}


// ************************************************************************* //
