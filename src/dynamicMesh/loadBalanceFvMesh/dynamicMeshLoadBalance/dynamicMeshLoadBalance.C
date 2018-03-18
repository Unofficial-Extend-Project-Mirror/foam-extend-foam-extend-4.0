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

#include "dynamicMeshLoadBalance.H"
#include "volFields.H"

#include "domainDecomposition.H"
#include "fvFieldDecomposer.H"

#include "processorMeshesReconstructor.H"
#include "fvFieldReconstructor.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// namespace Foam
// {
//     defineTypeNameAndDebug(dynamicMeshLoadBalance, 0);

//     addToRunTimeSelectionTable
//     (
//         topoChangerFvMesh,
//         dynamicFvMesh,
//         IOobject
//     );
// }



// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicMeshLoadBalance::dynamicMeshLoadBalance(dynamicFvMesh& mesh)
:
    mesh_(mesh),
    dict_
    (
        IOdictionary
        (
            IOobject
            (
                "decomposeParDict",
                mesh_.time().system(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        )
    ),
    imbalanceTrigger_(readScalar(dict_.lookup("imbalanceTrigger")))
{
    // Check imbalance trigger
    if (imbalanceTrigger_ < SMALL || imbalanceTrigger_ > 1)
    {
        WarningIn
        (
            "dynamicMeshLoadBalance::"
            "dynamicMeshLoadBalance(dynamicFvMesh& mesh)"
        )   << "Invalid imbalance trigger " << imbalanceTrigger_
            << " Should be between 0 and 1.  Resetting to 0.8"
            << endl;

        imbalanceTrigger_ = 0.8;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicMeshLoadBalance::~dynamicMeshLoadBalance()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mapPolyMesh>
Foam::dynamicMeshLoadBalance::updateMesh() const
{
    Info<< "Hello from dynamicMeshLoadBalance::updateMesh()" << endl;
    
    // Create a parallel decomposition
    domainDecomposition meshDecomp
    (
        mesh_,
        dict_
    );

    // Decompose the mesh based on processor data
    meshDecomp.decomposeMesh(false);

    // Analyse mesh decomposition to see how many cells are passed to which processor
    labelListList migratedCells(meshDecomp.nProcs());

    // Fill local list with number of local cells to be sent to each processor
    // First index = my processor.  Second index = processor to send to
    labelList& curMigratedCells = migratedCells[Pstream::myProcNo()];

    curMigratedCells.setSize(Pstream::myProcNo(), 0);

    const labelList& cellToProc = meshDecomp.cellToProc();

    forAll (cellToProc, cellI)
    {
        curMigratedCells[cellToProc[cellI]]++;
    }

    // Now that each processor has filled in its own part, combine the data
    Pstream::gatherList(migratedCells);
    Pstream::scatterList(migratedCells);

    // Reading through second index now tells how many cells will arrive
    // from which processor


    // Split the mesh and fields over processors and send
    for (label procI = 0; procI < meshDecomp.nProcs(); procI++)
    {
        if (procI != Pstream::myProcNo())
        {
            Info<< "Processor " << procI << ": send mesh and fields"
                << endl;

            // Check if there is a mesh to send
            if (curMigratedCells[procI] > 0)
            {
                // Send mesh and field to processor procI

                // Get mesh from decomposition
                autoPtr<fvMesh> procMeshPtr = meshDecomp.processorMesh
                (
                    procI,
                    mesh_.time(),
                    "processorPart" + name(procI)
                );
                fvMesh& procMesh = procMeshPtr();

                // Create a field decomposer
                fvFieldDecomposer fieldDecomposer
                (
                    mesh_,
                    procMesh,
                    meshDecomp.procFaceAddressing()[procI],
                    meshDecomp.procCellAddressing()[procI],
                    meshDecomp.procBoundaryAddressing()[procI]
                );
/*
                // Simple test only: rebalance volScalarFields
                HashTable<const volScalarField*> volScalarFields =
                    mesh_.thisDb().lookupClass<volScalarField>();

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

                // Note: communication can be optimised.  HJ, 27/Feb/2018
                OPstream toProc
                (
                    Pstream::nonBlocking,
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
                    Pstream::nonBlocking,
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
                            "processorPart" + name(procI),
                            mesh_.time().timeName(),
                            mesh_.time(),
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        fromProc,
                        false             // Do not sync
                    )
                );

                // Receive the fields
            }
        }
        else
        {
            // Insert own mesh
        }
    }

    // Create the reconstructed mesh
    autoPtr<fvMesh> reconstructedMeshPtr =
        meshRecon.reconstructMesh(mesh_.time());
    const fvMesh& reconstructedMesh = reconstructedMeshPtr();

    // Create field reconstructor
    fvFieldReconstructor fieldReconstructor
    (
        mesh_,
        procMeshes,
        meshRecon.faceProcAddressing(),
        meshRecon.cellProcAddressing(),
        meshRecon.boundaryProcAddressing()
    );

    // Create mapPolyMesh

    // DUMMY!!!
    autoPtr<mapPolyMesh> tmpm(new mapPolyMesh(mesh_));

    return tmpm;
}


// ************************************************************************* //
