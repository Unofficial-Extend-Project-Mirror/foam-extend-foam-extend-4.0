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
            
    Info<< "Hello from loadBalanceFvMesh::update()" << endl;

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

    // Analyse mesh decomposition to see how many cells are passed to which processor
    labelListList migratedCells(meshDecomp.nProcs());

    // Fill local list with number of local cells to be sent to each processor
    // First index = my processor.  Second index = processor to send to
    labelList& curMigratedCells = migratedCells[Pstream::myProcNo()];

    curMigratedCells.setSize(Pstream::nProcs(), 0);

    const labelList& cellToProc = meshDecomp.cellToProc();

    forAll (cellToProc, cellI)
    {
        curMigratedCells[cellToProc[cellI]]++;
    }

    // Now that each processor has filled in its own part, combine the data
    Pstream::gatherList(migratedCells);
    Pstream::scatterList(migratedCells);
    Pout<< "migratedCells: " << migratedCells << endl;
    // Reading through second index now tells how many cells will arrive
    // from which processor


    // Split the mesh and fields over processors and send
    for (label procI = 0; procI < meshDecomp.nProcs(); procI++)
    {
        if (procI != Pstream::myProcNo())
        {
            Pout<< "Send mesh and fields to processor " << procI << endl;

            // Check if there is a mesh to send
            if (curMigratedCells[procI] > 0)
            {
                // Send mesh and field to processor procI

                // Get mesh from decomposition
                autoPtr<fvMesh> procMeshPtr = meshDecomp.processorMesh
                (
                    procI,
                    time(),
                    "processorPart" + Foam::name(procI)
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
/*
                // Simple test only: rebalance volScalarFields
                HashTable<const volScalarField*> volScalarFields =
                    thisDb().lookupClass<volScalarField>();

                for
                (
                    HashTable<const volScalarField*>::const_iterator iter =
                        volScalarFields.begin();
                    iter != volScalarFields.end();
                    ++iter
                )
                {
                    fieldDecomposer.decomposeField(*(iter()));
                }
*/

                OPstream toProc
                (
                    Pstream::blocking,
                    procI
                );

                // Send the mesh and fields to target processor
                toProc << procMesh;
            }
        }
    }


    // Collect pieces of mesh and fields from other processors

    // Create a list of processor meshes
    // Create the reconstructor
    processorMeshesReconstructor meshRecon("reconstructed");

    PtrList<fvMesh>& procMeshes = meshRecon.meshes();
    procMeshes.setSize(meshDecomp.nProcs());

    List<PtrList<volScalarField*> > receiveVolScalarFields
    (
        meshDecomp.nProcs()
    );

    for (label procI = 0; procI < meshDecomp.nProcs(); procI++)
    {
        if (procI != Pstream::myProcNo())
        {
            // Check if there is a mesh to send
            if (migratedCells[procI][Pstream::myProcNo()] > 0)
            {
                Info<< "Processor " << procI << ": receive mesh and fields"
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
                Pout<< "Received " << procI << " mesh stats: "
                    << " polyPatches: "
                    << procMeshes[procI].boundaryMesh().size()
                    << " patches: "
                    << procMeshes[procI].boundary().size()
                    << endl;
                // Receive the fields
            }
        }
        else
        {
            // Insert own mesh
            procMeshes.set(procI, this);
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

    // Create field reconstructor
    // fvFieldReconstructor fieldReconstructor
    // (
    //     *this,
    //     procMeshes,
    //     meshRecon.faceProcAddressing(),
    //     meshRecon.cellProcAddressing(),
    //     meshRecon.boundaryProcAddressing()
    // );


    // Note: Local mesh and fields need to be removed from the PtrList
    // before the destruction and from autoPtr.
    // Otherwise, they will be deleted.  HJ, 5/Mar/2018
    autoPtr<fvMesh> curMesh = procMeshes.set(Pstream::myProcNo(), NULL);
    curMesh.ptr();
    
    return true;
}


// ************************************************************************* //
