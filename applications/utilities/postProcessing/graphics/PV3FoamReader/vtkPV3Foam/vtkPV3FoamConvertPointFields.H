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

#ifndef vtkPV3FoamConvertPointFields_H
#define vtkPV3FoamConvertPointFields_H

// Foam includes
#include "interpolatePointToCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPV3Foam::convertPointFields
(
    const fvMesh& mesh,
    const IOobjectList& objects,
    vtkDataArraySelection* fieldSelection,
    vtkMultiBlockDataSet* output
)
{
    IOobjectList fieldObjects
    (
        objects.lookupClass
        (
            GeometricField<Type, pointPatchField, pointMesh>::typeName
        )
    );

    label nSelectedFields = fieldSelection->GetNumberOfArrays();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    for (label i=0; i<nSelectedFields; i++)
    {
        if (fieldSelection->GetArraySetting(i))
        {
            word fieldName = fieldSelection->GetArrayName(i);

            if (fieldObjects.found(fieldName))
            {

                if (debug)
                {
                    Info<< "Foam::vtkPV3Foam::convertPointFields : "
                        << fieldName << endl;
                }

                pointMesh pMesh(mesh);

                GeometricField<Type, pointPatchField, pointMesh> ptf
                (
                    IOobject
                    (
                        fieldName,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    pMesh
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
                        convertPointField
                        (
                            ptf,
                            GeometricField<Type, fvPatchField, volMesh>::null(),
                            output,
                            selectInfoVolume_,
                            selectedRegionDatasetIds_[regionId]
                        );
                    }
                }


                // Convert patches
                label regionId = selectInfoPatches_.start();
                forAll (patches, patchI)
                {
                    if (patches[patchI].size())
                    {
                        if (selectedRegions_[regionId])
                        {
                            convertPatchPointField
                            (
                                fieldName,
                                ptf.boundaryField()[patchI]
                                   .patchInternalField()(),
                                output,
                                selectInfoPatches_,
                                selectedRegionDatasetIds_[regionId]
                            );
                        }
                        regionId++;
                    }
                }
            }
        }
    }
}


template<class Type>
void Foam::vtkPV3Foam::convertPointField
(
    const GeometricField<Type, pointPatchField, pointMesh>& ptf,
    const GeometricField<Type, fvPatchField, volMesh>& tf,
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector,
    const label datasetNo
)
{
    const label nComp = pTraits<Type>::nComponents;

    vtkUnstructuredGrid* internalMesh = vtkUnstructuredGrid::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    );

    vtkFloatArray *pointData = vtkFloatArray::New();

    // Recalculate size of data field to include inactive points (supporting
    // inactive faces).  Data will be padded for correct size
    // HJ, 28/Mar/2009
    pointData->SetNumberOfTuples
    (
        ptf.size()
      + addPointCellLabels_.size()
      + tf.mesh().allPoints().size() - tf.mesh().nPoints()
    );
    pointData->SetNumberOfComponents(nComp);
    pointData->Allocate(nComp*ptf.size());
    pointData->SetName(tf.name().c_str());

    float vec[nComp];

    forAll(ptf, i)
    {
        for (direction d=0; d<nComp; d++)
        {
            vec[d] = component(ptf[i], d);
        }

        pointData->InsertTuple(i, vec);
    }

    label i = ptf.size();

    if (&tf != &GeometricField<Type, fvPatchField, volMesh>::null())
    {
        forAll(addPointCellLabels_, api)
        {
            Type t = tf[addPointCellLabels_[api]];

            for (direction d=0; d<nComp; d++)
            {
                vec[d] = component(t, d);
            }

            pointData->InsertTuple(i++, vec);
        }
    }
    else
    {
        forAll(addPointCellLabels_, api)
        {
            Type t = interpolatePointToCell(ptf, addPointCellLabels_[api]);

            for (direction d=0; d<nComp; d++)
            {
                vec[d] = component(t, d);
            }

            pointData->InsertTuple(i++, vec);
        }
    }

    // Pad for inactive points
    // HJ, 28/Mar/2009
    for (direction d = 0; d < nComp; d++)
    {
        vec[d] = 0;
    }

    for
    (
        label padI = tf.mesh().nPoints();
        padI < tf.mesh().allPoints().size();
        padI++
    )
    {
        pointData->InsertTuple(i++, vec);
    }

    internalMesh->GetPointData()->AddArray(pointData);
    pointData->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
