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

InClass
    vtkPV3Foam

\*---------------------------------------------------------------------------*/

#ifndef vtkPV3FoamConvertVolFields_H
#define vtkPV3FoamConvertVolFields_H

// Foam includes
#include "emptyFvPatchField.H"
#include "wallPolyPatch.H"
#include "faceSet.H"
#include "vtkPV3FoamConvertPatchFaceField.H"
#include "vtkPV3FoamConvertPatchPointField.H"
#include "vtkPV3FoamConvertFaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPV3Foam::convertVolFields
(
    const fvMesh& mesh,
    const volPointInterpolation& pInterp,
    const PtrList<PrimitivePatchInterpolation<primitivePatch> >& ppInterpList,
    const IOobjectList& objects,
    vtkDataArraySelection *fieldSelection,
    vtkMultiBlockDataSet* output
)
{
    IOobjectList fieldObjects
    (
        objects.lookupClass
        (
            GeometricField<Type, fvPatchField, volMesh>::typeName
        )
    );

    label nFields = fieldSelection->GetNumberOfArrays();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    vtkDataArraySelection* arraySelector = reader_->GetRegionSelection();

    for (label i=0; i<nFields; i++)
    {
        const word fieldName = fieldSelection->GetArrayName(i);

        if
        (
            !fieldSelection->GetArraySetting(i)
         || !fieldObjects.found(fieldName))
        {
            continue;
        }

        if (debug)
        {
            Info<< "converting Foam volume field: " << fieldName
                << endl;
        }

        GeometricField<Type, fvPatchField, volMesh> tf
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        tmp<GeometricField<Type, pointPatchField, pointMesh> > tptf
        (
            pInterp.interpolate(tf)
        );

        // Convert internal mesh
        for
        (
            int regionId = selectInfoVolume_.start();
            regionId < selectInfoVolume_.end();
            ++regionId
        )
        {
            if (selectedRegions_[regionId])
            {
                convertVolField
                (
                    tf,
                    output,
                    selectInfoVolume_,
                    selectedRegionDatasetIds_[regionId],
                    superCells_
                );
                convertPointField
                (
                    tptf(),
                    tf,
                    output,
                    selectInfoVolume_,
                    selectedRegionDatasetIds_[regionId]
                );
            }
        }


        // Convert patches
        for
        (
            int regionId = selectInfoPatches_.start();
            regionId < selectInfoPatches_.end();
            ++regionId
        )
        {
            if (!selectedRegions_[regionId])
            {
                continue;
            }

            word selectName = getFirstWord
            (
                arraySelector->GetArrayName(regionId)
            );

            const label patchId = patches.findPatchID(selectName);

            const fvPatchField<Type>& ptf
            (
                tf.boundaryField()[patchId]
            );

            if
            (
                isType<emptyFvPatchField<Type> >(ptf)
             ||
                (
                    typeid(patches[patchId]) == typeid(wallPolyPatch)
                 && reader_->GetExtrapolateWalls()
                )
            )
            {
                fvPatch p(ptf.patch().patch(), tf.mesh().boundary());
                tmp<Field<Type> > tpptf
                (
                    fvPatchField<Type>(p, tf).patchInternalField()
                );

                convertPatchFaceField
                (
                    tf.name(),
                    tpptf(),
                    output,
                    selectInfoPatches_,
                    selectedRegionDatasetIds_[regionId]
                );

                convertPatchPointField
                (
                    tf.name(),
                    ppInterpList[patchId].faceToPointInterpolate(tpptf)(),
                    output,
                    selectInfoPatches_,
                    selectedRegionDatasetIds_[regionId]
                );
            }
            else
            {
                convertPatchFaceField
                (
                    tf.name(),
                    ptf,
                    output,
                    selectInfoPatches_,
                    selectedRegionDatasetIds_[regionId]
                );
                convertPatchPointField
                (
                    tf.name(),
                    ppInterpList[patchId].faceToPointInterpolate(ptf)(),
                    output,
                    selectInfoPatches_,
                    selectedRegionDatasetIds_[regionId]
                );
            }
        }

        // Convert cell zones
        for
        (
            int regionId = selectInfoCellZones_.start();
            regionId < selectInfoCellZones_.end();
            ++regionId
        )
        {
            if (!selectedRegions_[regionId])
            {
                continue;
            }

            if (debug)
            {

                word selectName = getFirstWord
                (
                    arraySelector->GetArrayName(regionId)
                );

                Info<< "wish to convert cellzone: " << selectName
                    << " regionId: " << regionId
                    << " volume field: " << fieldName
                    << endl;
            }

            const label datasetId =
                selectedRegionDatasetIds_[regionId];

            convertVolField
            (
                tf,
                output, selectInfoCellZones_, datasetId,
                zoneSuperCells_[datasetId]
            );
        }

        // Convert cell sets
        for
        (
            int regionId = selectInfoCellSets_.start();
            regionId < selectInfoCellSets_.end();
            ++regionId
        )
        {
            if (!selectedRegions_[regionId])
            {
                continue;
            }

            if (debug)
            {
                word selectName = getFirstWord
                (
                    arraySelector->GetArrayName(regionId)
                );

                Info<< "wish to convert cellset: " << selectName
                    << " regionId: " << regionId
                    << " volume field: " << fieldName
                    << endl;
            }

            const label datasetId =
                selectedRegionDatasetIds_[regionId];

            convertVolField
            (
                tf,
                output, selectInfoCellSets_, datasetId,
                csetSuperCells_[datasetId]
            );
        }

        // Convert face zones
        for
        (
            int regionId = selectInfoFaceZones_.start();
            regionId < selectInfoFaceZones_.end();
            ++regionId
        )
        {
            if (!selectedRegions_[regionId])
            {
                continue;
            }

            if (debug)
            {
                word selectName = getFirstWord
                (
                    arraySelector->GetArrayName(regionId)
                );

                Info<< "wish to convert facezone: " << selectName
                    << " regionId: " << regionId
                    << " volume field: " << fieldName
                    << endl;
            }

            const faceZoneMesh& fzMesh = mesh.faceZones();
            const label zoneI = regionId - selectInfoFaceZones_.start();

            const label datasetId =
                selectedRegionDatasetIds_[regionId];

            convertFaceField
            (
                tf,
                output, selectInfoFaceZones_, datasetId,
                mesh,
                fzMesh[zoneI]
            );
        }

        // Convert face sets
        for
        (
            int regionId = selectInfoFaceSets_.start();
            regionId < selectInfoFaceSets_.end();
            ++regionId
        )
        {
            if (!selectedRegions_[regionId])
            {
                continue;
            }

            word selectName = getFirstWord
            (
                arraySelector->GetArrayName(regionId)
            );

            if (debug)
            {
                Info<< "wish to convert faceset: " << selectName
                    << " regionId: " << regionId
                    << " volume field: " << fieldName
                    << endl;
            }

            const faceSet fSet(mesh, selectName);

            const label datasetId =
                selectedRegionDatasetIds_[regionId];

            convertFaceField
            (
                tf,
                output, selectInfoFaceSets_, datasetId,
                mesh,
                fSet
            );
        }
    }
}


template<class Type>
void Foam::vtkPV3Foam::convertVolField
(
    const GeometricField<Type, fvPatchField, volMesh>& tf,
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector,
    const label datasetNo,
    labelList& superCells
)
{
    const label nComp = pTraits<Type>::nComponents;

    vtkUnstructuredGrid* vtkmesh = vtkUnstructuredGrid::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    );

    vtkFloatArray* celldata = vtkFloatArray::New();
    celldata->SetNumberOfTuples(superCells.size());
    celldata->SetNumberOfComponents(nComp);
    celldata->Allocate(nComp*superCells.size());
    celldata->SetName(tf.name().c_str());

    if (debug)
    {
        Info<< "converting vol<Type>Field: " << tf.name() << nl
            << "field size = " << tf.size() << nl
            << "nTuples = " << superCells.size() << nl
            << "nComp  = " << nComp << endl;
    }

    float vec[nComp];

    forAll(superCells, scI)
    {
        const Type& t = tf[superCells[scI]];
        for (direction d=0; d<nComp; d++)
        {
            vec[d] = component(t, d);
        }

        celldata->InsertTuple(scI, vec);
    }

    vtkmesh->GetCellData()->AddArray(celldata);
    celldata->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
