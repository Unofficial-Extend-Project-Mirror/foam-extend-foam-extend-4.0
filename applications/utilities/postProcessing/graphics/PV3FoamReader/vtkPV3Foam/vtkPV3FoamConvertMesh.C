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

Description

\*---------------------------------------------------------------------------*/

#include "vtkPV3Foam.H"

// Foam includes
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "fvMeshSubset.H"
#include "vtkPV3FoamReader.h"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPV3Foam::convertMeshVolume
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertMeshVolume" << endl;
        printMemory();
    }

    const selectionInfo& selector = selectInfoVolume_;
    const fvMesh& mesh = *meshPtr_;

    // Create the internal mesh and add as dataset 0
    for
    (
        int regionId = selector.start();
        regionId < selector.end();
        ++regionId
    )
    {
        if (!selectedRegions_[regionId])
        {
            continue;
        }

        // word selectName = getFirstWord
        // (
        //     arraySelection->GetArrayName(regionId)
        // );

        if (debug)
        {
            Info<< "Creating VTK internalMesh" << endl;
        }

        const label datasetId = 0;

        vtkUnstructuredGrid* vtkmesh = vtkUnstructuredGrid::New();
        addVolumeMesh(mesh, vtkmesh, superCells_);

        AddToBlock(output, selector, datasetId, vtkmesh, "internalMesh");
        selectedRegionDatasetIds_[regionId] = datasetId;
        vtkmesh->Delete();
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertMeshVolume" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::convertMeshLagrangian
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertMeshLagrangian" << endl;
        printMemory();
    }

    const selectionInfo& selector = selectInfoLagrangian_;
    const fvMesh& mesh = *meshPtr_;

    // Create the Lagrangian mesh and add as dataset 0
    for
    (
        int regionId = selector.start();
        regionId < selector.end();
        ++regionId
    )
    {
        if (!selectedRegions_[regionId])
        {
            continue;
        }

        if (debug)
        {
            Info<< "Creating VTK Lagrangian mesh" << endl;
        }

        const label datasetId = 0;

        vtkPolyData* vtkmesh = vtkPolyData::New();
        addLagrangianMesh(mesh, vtkmesh);

        AddToBlock(output, selector, datasetId, vtkmesh, cloudName_);
        selectedRegionDatasetIds_[regionId] = datasetId;
        vtkmesh->Delete();
    }
    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertMeshLagrangian" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::convertMeshPatches
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertMeshPatches" << endl;
        printMemory();
    }

    const selectionInfo& selector = selectInfoPatches_;
    vtkDataArraySelection* arraySelection = reader_->GetRegionSelection();

    if (selector.size())
    {
        const fvMesh& mesh = *meshPtr_;
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        // Create the patches and add as dataset ...
        for
        (
            int regionId = selector.start();
            regionId < selector.end();
            ++regionId
        )
        {
            if (!selectedRegions_[regionId])
            {
                continue;
            }

            word selectName = getFirstWord
            (
                arraySelection->GetArrayName(regionId)
            );

            const label patchId = patches.findPatchID(selectName);

            if (debug)
            {
                Info<< "Creating VTK mesh for patch: " << selectName
                    << " region index: " << regionId << endl;
            }

            const label datasetId = GetNumberOfDataSets(output, selector);

            vtkPolyData* vtkmesh = vtkPolyData::New();
            addPatchMesh
            (
                patches[patchId],
                vtkmesh
            );

            AddToBlock
            (
                output, selector, datasetId, vtkmesh,
                selectName + ":patch"
            );
            selectedRegionDatasetIds_[regionId] = datasetId;
            vtkmesh->Delete();
        }
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertMeshPatches" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::convertMeshCellZones
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertMeshCellZones" << endl;
        printMemory();
    }

    const selectionInfo& selector = selectInfoCellZones_;
    const fvMesh& mesh = *meshPtr_;

    // Create the cell zone(s) and add as DataSet(CELLZONE, 0..n)
    if (selector.size())
    {
        const cellZoneMesh& czMesh = mesh.cellZones();

        // use the zoneId directly instead of the name
        for (int zoneI=0; zoneI < selector.size(); ++zoneI)
        {
            const int regionId = selector.start() + zoneI;

            if (!selectedRegions_[regionId])
            {
                continue;
            }

            if (debug)
            {
                Info<< "Creating VTK mesh for cellZone: "
                    << zoneI << endl;
            }

            fvMeshSubset subsetter
            (
                IOobject
                (
                    "set",
                    mesh.time().constant(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            subsetter.setLargeCellSubset(labelHashSet(czMesh[zoneI]));

            const label datasetId = GetNumberOfDataSets(output, selector);

            vtkUnstructuredGrid* vtkmesh = vtkUnstructuredGrid::New();

            addVolumeMesh
            (
                subsetter.subMesh(),
                vtkmesh,
                zoneSuperCells_[datasetId]
            );

            // renumber - superCells must contain global cell ids
            inplaceRenumber
            (
                subsetter.cellMap(),
                zoneSuperCells_[datasetId]
            );

            AddToBlock
            (
                output, selector, datasetId, vtkmesh,
                czMesh.names()[zoneI] + ":cellZone"
            );
            selectedRegionDatasetIds_[regionId] = datasetId;
            vtkmesh->Delete();
        }
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertMeshCellZones" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::convertMeshCellSets
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertMeshCellSets" << endl;
        printMemory();
    }

    const selectionInfo& selector = selectInfoCellSets_;
    vtkDataArraySelection* arraySelection = reader_->GetRegionSelection();

    // Create the cell sets and add as dataset
    if (selector.size())
    {
        const fvMesh& mesh = *meshPtr_;

        for
        (
            int regionId = selector.start();
            regionId < selector.end();
            ++regionId)
        {
            if (!selectedRegions_[regionId])
            {
                continue;
            }

            word selectName = getFirstWord
            (
                arraySelection->GetArrayName(regionId)
            );

            if (debug)
            {
                Info<< "Creating VTK mesh for cellSet: " << selectName
                    << " region index: " << regionId << endl;
            }

            const cellSet cSet(mesh, selectName);
            fvMeshSubset subsetter
            (
                IOobject
                (
                    "set",
                    mesh.time().constant(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );
            subsetter.setLargeCellSubset(cSet);

            const label datasetId = GetNumberOfDataSets(output, selector);

            vtkUnstructuredGrid* vtkmesh = vtkUnstructuredGrid::New();

            addVolumeMesh
            (
                subsetter.subMesh(),
                vtkmesh,
                csetSuperCells_[datasetId]
            );

            // renumber - superCells must contain global cell ids
            inplaceRenumber
            (
                subsetter.cellMap(),
                csetSuperCells_[datasetId]
            );

            AddToBlock
            (
                output, selector, datasetId, vtkmesh,
                selectName + ":cellSet"
            );
            selectedRegionDatasetIds_[regionId] = datasetId;
            vtkmesh->Delete();
        }
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertMeshCellSets" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::convertMeshFaceZones
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertMeshFaceZones" << endl;
        printMemory();
    }

    const selectionInfo& selector = selectInfoFaceZones_;
    const fvMesh& mesh = *meshPtr_;

    // Create the cell zone(s) and add as datasets
    if (selector.size())
    {
        const faceZoneMesh& fzMesh = mesh.faceZones();

        // use the zoneId directly instead of the name
        for (int zoneI=0; zoneI < selector.size(); ++zoneI)
        {
            const int regionId = selector.start() + zoneI;

            if (!selectedRegions_[regionId])
            {
                continue;
            }

            if (debug)
            {
                Info<< "Creating VTK mesh for faceZone: "
                    << zoneI << endl;
            }

            const label datasetId = GetNumberOfDataSets(output, selector);

            vtkPolyData* vtkmesh = vtkPolyData::New();

            addFaceZoneMesh
            (
                mesh,
                fzMesh[zoneI],
                vtkmesh
            );

            AddToBlock
            (
                output, selector, datasetId, vtkmesh,
                fzMesh.names()[zoneI] + ":faceZone"
            );
            selectedRegionDatasetIds_[regionId] = datasetId;
            vtkmesh->Delete();
        }
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertMeshFaceZones" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::convertMeshFaceSets
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertMeshFaceSets" << endl;
        printMemory();
    }

    const selectionInfo& selector = selectInfoFaceSets_;
    vtkDataArraySelection* arraySelection = reader_->GetRegionSelection();

    // Create the face sets and add as dataset
    if (selector.size())
    {
        const fvMesh& mesh = *meshPtr_;

        for
        (
            int regionId = selector.start();
            regionId < selector.end();
            ++regionId
        )
        {
            if (!selectedRegions_[regionId])
            {
                continue;
            }

            word selectName = getFirstWord
            (
                arraySelection->GetArrayName(regionId)
            );

            if (debug)
            {
                Info<< "Creating VTK mesh for faceSet: " << selectName
                    << " region index: " << regionId << endl;
            }

            const faceSet fSet(mesh, selectName);

            const label datasetId = GetNumberOfDataSets(output, selector);

            vtkPolyData* vtkmesh = vtkPolyData::New();
            addFaceSetMesh
            (
                mesh,
                fSet,
                vtkmesh
            );

            AddToBlock
            (
                output, selector, datasetId, vtkmesh,
                selectName + ":faceSet"
            );
            selectedRegionDatasetIds_[regionId] = datasetId;
            vtkmesh->Delete();
        }
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertMeshFaceSets" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::convertMeshPointZones
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertMeshPointZones" << endl;
        printMemory();
    }

    const selectionInfo& selector = selectInfoPointZones_;
    const fvMesh& mesh = *meshPtr_;

    // Create the point sets and add as dataset
    if (selector.size())
    {
        const pointZoneMesh& pzMesh = mesh.pointZones();

        // use the zoneId directly instead of the name
        for (int zoneI=0; zoneI < selector.size(); ++zoneI)
        {
            const int regionId = selector.start() + zoneI;

            if (!selectedRegions_[regionId])
            {
                continue;
            }

            if (debug)
            {
                Info<< "Creating VTK mesh for pointZone: "
                    << zoneI << endl;
            }

            const label datasetId = GetNumberOfDataSets(output, selector);

            vtkPolyData* vtkmesh = vtkPolyData::New();

            addPointZoneMesh
            (
                mesh,
                pzMesh[zoneI],
                vtkmesh
            );

            AddToBlock
            (
                output, selector, datasetId, vtkmesh,
                pzMesh.names()[zoneI] + ":pointZone"
            );
            selectedRegionDatasetIds_[regionId] = datasetId;
            vtkmesh->Delete();
        }
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertMeshPointZones" << endl;
        printMemory();
    }
}



void Foam::vtkPV3Foam::convertMeshPointSets
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertMeshPointSets" << endl;
        printMemory();
    }

    const selectionInfo& selector = selectInfoPointSets_;
    vtkDataArraySelection* arraySelection = reader_->GetRegionSelection();

    // Create the point sets and add as dataset
    if (selector.size())
    {
        const fvMesh& mesh = *meshPtr_;

        for
        (
            int regionId = selector.start();
            regionId < selector.end();
            ++regionId
        )
        {
            if (!selectedRegions_[regionId])
            {
                continue;
            }

            word selectName = getFirstWord
            (
                arraySelection->GetArrayName(regionId)
            );


            if (debug)
            {
                Info<< "Creating VTK mesh for pointSet: " << selectName
                    << " region index: " << regionId << endl;
            }

            const pointSet pSet(mesh, selectName);

            const label datasetId = GetNumberOfDataSets(output, selector);

            vtkPolyData* vtkmesh = vtkPolyData::New();
            addPointSetMesh
            (
                mesh,
                pSet,
                vtkmesh
            );
            AddToBlock
            (
                output, selector, datasetId, vtkmesh,
                selectName + ":pointSet"
            );
            selectedRegionDatasetIds_[regionId] = datasetId;
            vtkmesh->Delete();
        }
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertMeshPointSets" << endl;
        printMemory();
    }
}

// ************************************************************************* //
