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

#ifndef vtkPV3FoamConvertLagrangianFields_H
#define vtkPV3FoamConvertLagrangianFields_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPV3Foam::convertLagrangianFields
(
    const fvMesh& mesh,
    const IOobjectList& objects,
    vtkDataArraySelection* fieldSelection,
    vtkMultiBlockDataSet* output
)
{
    const selectionInfo& selector = selectInfoLagrangian_;
    const label datasetId = 0;

    if (!selector.size())
    {
        return;
    }

    IOobjectList fieldObjects
    (
        objects.lookupClass
        (
            IOField<Type>::typeName
        )
    );

    label nFields = fieldSelection->GetNumberOfArrays();

    for (label i=0; i<nFields; i++)
    {
        if (fieldSelection->GetArraySetting(i))
        {
            const word fieldName = fieldSelection->GetArrayName(i);

            if (fieldObjects.found(fieldName))
            {
                IOField<Type> iotf
                (
                    IOobject
                    (
                        fieldName,
                        mesh.time().timeName(),
                        "lagrangian"/cloudName_,
                        mesh,
                        IOobject::MUST_READ
                    )
                );

                convertLagrangianField(iotf, output, selector, datasetId);
            }
        }
    }
}


template<class Type>
void Foam::vtkPV3Foam::convertLagrangianField
(
    const IOField<Type>& tf,
    vtkMultiBlockDataSet* output,
    const selectionInfo& selector,
    const label datasetNo
)
{
    const label nComp = pTraits<Type>::nComponents;

    vtkPolyData* vtkmesh = vtkPolyData::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    );

    vtkFloatArray *pointData = vtkFloatArray::New();
    pointData->SetNumberOfTuples(tf.size());
    pointData->SetNumberOfComponents(nComp);
    pointData->Allocate(nComp*tf.size());
    pointData->SetName(tf.name().c_str());

    if (debug)
    {
        Info<< "converting Lagrangian <Type>Field: " << tf.name() << nl
            << "tf.size() = " << tf.size() << nl
            << "nComp  = " << nComp << endl;
    }

    float vec[nComp];

    forAll(tf, i)
    {
        for (direction d=0; d<nComp; d++)
        {
            vec[d] = component(tf[i], d);
        }

        pointData->InsertTuple(i, vec);
    }

    vtkmesh->GetPointData()->AddArray(pointData);
    pointData->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
