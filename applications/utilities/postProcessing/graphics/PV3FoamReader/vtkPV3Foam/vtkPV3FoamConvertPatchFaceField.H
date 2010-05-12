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

#ifndef vtkPV3FoamConvertPatchFaceField_H
#define vtkPV3FoamConvertPatchFaceField_H

// VTK includes
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPV3Foam::convertPatchFaceField
(
    const word& name,
    const Field<Type>& ptf,
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

    vtkFloatArray *cellData = vtkFloatArray::New();
    cellData->SetNumberOfTuples(ptf.size());
    cellData->SetNumberOfComponents(nComp);
    cellData->Allocate(nComp*ptf.size());
    cellData->SetName(name.c_str());

    float vec[nComp];

    forAll(ptf, i)
    {
        const Type& t = ptf[i];
        for (direction d=0; d<nComp; d++)
        {
            vec[d] = component(t, d);
        }

        cellData->InsertTuple(i, vec);
    }

    vtkmesh->GetCellData()->AddArray(cellData);
    cellData->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
