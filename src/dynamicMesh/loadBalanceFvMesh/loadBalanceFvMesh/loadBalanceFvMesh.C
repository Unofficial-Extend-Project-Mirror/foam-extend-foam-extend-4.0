/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "loadBalanceFvMesh.H"
#include "domainDecomposition.H"
#include "fvFieldDecomposer.H"
#include "processorMeshesReconstructor.H"
#include "fvFieldReconstructor.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"

#include "GeoMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(loadBalanceFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, loadBalanceFvMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::loadBalanceFvMesh::checkLoadBalance
(
    const scalarField& weights
) const
{
    if (Pstream::parRun())
    {
        // Calculate local and global load
        scalar localLoad = sum(weights);

        scalar globalLoad = localLoad;

        reduce(globalLoad, sumOp<scalar>());

        globalLoad /= Pstream::nProcs();

        // Calculate imbalance as min of localLoad/globalLoad
        scalar imbalance = mag(1 - localLoad/globalLoad);

        reduce(imbalance, minOp<scalar>());

        Info<< "Global imbalance: " << imbalance << endl;

        if (imbalance < imbalanceTrigger_)
        {
            return true;
        }
    }

    // Serial run or low imbalance: no balancing possible
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::loadBalanceFvMesh::loadBalanceFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    imbalanceTrigger_
    (
        readScalar(dynamicMeshCoeffs_.lookup("imbalanceTrigger"))
    )
{
    // Check imbalance trigger
    if (imbalanceTrigger_ < SMALL || imbalanceTrigger_ > 1)
    {
        WarningIn
        (
            "loadBalanceFvMesh::"
            "loadBalanceFvMesh(const IOobject& io)"
        )   << "Invalid imbalance trigger " << imbalanceTrigger_
            << " Should be between 0 and 1.  Resetting to 0.8"
            << endl;

        imbalanceTrigger_ = 0.8;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::loadBalanceFvMesh::~loadBalanceFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::loadBalanceFvMesh::update()
{
    if (!Pstream::parRun())
    {
        InfoIn("bool loadBalanceFvMesh::update()")
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
        dynamicMeshCoeffs_
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

    // Write as volScalarField for post-processing
    volScalarField cellDist
    (
        IOobject
        (
            "cellDist",
            time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        dimensionedScalar("cellDist", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    );

    forAll(cellToProc, celli)
    {
       cellDist[celli] = cellToProc[celli];
    }

    cellDist.write();


    forAll (cellToProc, cellI)
    {
        curMigratedCells[cellToProc[cellI]]++;
    }

    // Now that each processor has filled in its own part, combine the data
    Pstream::gatherList(migratedCells);
    Pstream::scatterList(migratedCells);

    // Reading through second index now tells how many cells will arrive
    // from which processor

    // Find out which processor faces will become local internal faces
    // by comparing decomposition index from the other side


    // Split the mesh and fields over processors and send

    // Prepare receiving side

    // Create the reconstructor
    processorMeshesReconstructor meshRecon("reconstructed");

    PtrList<fvMesh>& procMeshes = meshRecon.meshes();
    procMeshes.setSize(meshDecomp.nProcs());

    // Collect local fields for decomposition
    clearOut();
    
    // Collect fields for load balancing
    HashTable<const volScalarField*> volScalarFields =
        thisDb().lookupClass<volScalarField>();
    Pout<< "volScalarFields: " << volScalarFields.sortedToc() << endl;
    HashTable<const volVectorField*>volVectorFields =
        thisDb().lookupClass<volVectorField>();
    Pout<< "volVectorFields: " << volVectorFields.sortedToc() << endl;
    HashTable<const surfaceScalarField*> surfaceScalarFields =
        thisDb().lookupClass<surfaceScalarField>();
    Pout<< "surfaceScalarFields: " << surfaceScalarFields.sortedToc() << endl;
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

    for (label procI = 0; procI < meshDecomp.nProcs(); procI++)
    {
        // Check if there is a mesh to send
        if (curMigratedCells[procI] > 0)
        {
            // Send mesh and field to processor procI

            // Get mesh from decomposition
            autoPtr<fvMesh> procMeshPtr = meshDecomp.processorMesh
            (
                procI,
                time(),
                "processorPart" + Foam::name(procI),
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
                Pout<< "Send mesh and fields to processor " << procI << endl;

                OPstream toProc
                (
                    Pstream::blocking,
                    procI
                );

                // Send the mesh and fields to target processor
                toProc << procMesh << nl;

                // Send fields
                sendFields(volScalarFields, fieldDecomposer, toProc);
                sendFields(volVectorFields, fieldDecomposer, toProc);
                sendFields(surfaceScalarFields, fieldDecomposer, toProc);
            }
            else
            {
                // My processor.  Grab mesh and fields
                // Note: procI == Pstream::myProcNo()
                // HJ, 26/Apr/2018

                // Set local mesh piece
                procMeshes.set
                (
                    procI,
                    procMeshPtr
                );

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
                Pout<< "Receive mesh and fields from " << procI
                    << endl;

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
                            "processorPart" + Foam::name(procI),
                            time().timeName(),
                            time(),
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        fromProc,
                        false             // Do not sync
                    )
                );

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
            }
        }
    }


    // Create the reconstructed mesh
    autoPtr<fvMesh> reconstructedMeshPtr =
        meshRecon.reconstructMesh(time());
    fvMesh& reconMesh = reconstructedMeshPtr();

    Pout<< "Reconstructed mesh stats: "
        << " nCells: " << reconMesh.nCells()
        << " polyPatches: "
        << reconMesh.boundaryMesh().size()
        << " patches: "
        << reconMesh.boundary().size()
        << endl;

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
            FatalErrorIn("bool loadBalanceFvMesh::update()")
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
                FatalErrorIn("bool loadBalanceFvMesh::update()")
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
                FatalErrorIn("bool loadBalanceFvMesh::update()")
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

    Pout<< "New mesh: points " << nPoints()
        << " faces: " << nFaces()
        << " internal: " << nInternalFaces()
        << " cells: " << nCells()
        << " patches: " << reconPatchSizes.size() << endl;

    forAll (procMeshes, procI)
    {
        Pout<< "procMesh " << procI
            << " points " << procMeshes[procI].nPoints()
            << " faces: " << procMeshes[procI].nFaces()
            << " internal: " << procMeshes[procI].nInternalFaces()
            << " cells: " << procMeshes[procI].nCells()
            << " patches: " << procMeshes[procI].boundary().size()
            << endl;
    }

    // Pout<< "RECON face: " << meshRecon.faceProcAddressing() << nl
    //     << " cell: " << meshRecon.cellProcAddressing() << nl
    //     << " boundary: " << meshRecon.boundaryProcAddressing() << nl
    //     << endl;


    // To Do: build a reconstructor from addressing data
    Pout<< nl << nl
        << "FINISHED CYCLE" << nl << nl << nl
        << endl;

    // Create a dummy mapping object to re-size all fields

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
        resetFvPatchFlag
    );

    rebuildFields
    (
        volVectorFields,
        fieldReconstructor,
        receivedVolVectorFields,
        resetFvPatchFlag
    );

    rebuildFields
    (
        surfaceScalarFields,
        fieldReconstructor,
        receivedSurfaceScalarFields,
        resetFvPatchFlag
    );

    return true;
}


// ************************************************************************* //
