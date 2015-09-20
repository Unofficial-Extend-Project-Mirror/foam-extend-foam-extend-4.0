/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

InClass
    vtkPV4Foam

\*---------------------------------------------------------------------------*/

#ifndef vtkPV4FoamFaceField_H
#define vtkPV4FoamFaceField_H

// VTK includes
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"

#include "vtkFoamTupleRemap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPV4Foam::convertFaceField
(
    const GeometricField<Type, fvPatchField, volMesh>& tf,
    vtkMultiBlockDataSet* output,
    const partInfo& selector,
    const label datasetNo,
    const fvMesh& mesh,
    const labelList& faceLabels
)
{
    const label nComp = pTraits<Type>::nComponents;
    const label nInternalFaces = mesh.nInternalFaces();
    const label nFaces = mesh.nFaces();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeigh = mesh.faceNeighbour();

    vtkFloatArray *cellData = vtkFloatArray::New();
    cellData->SetNumberOfTuples( faceLabels.size() );
    cellData->SetNumberOfComponents( nComp );
    cellData->Allocate( nComp*faceLabels.size() );
    cellData->SetName( tf.name().c_str() );

    if (debug)
    {
        Info<< "convert convertFaceField: "
            << tf.name()
            << " size = " << tf.size()
            << " nComp=" << nComp
            << " nTuples = " << faceLabels.size() <<  endl;
    }

    float vec[nComp];

    // for interior faces: average owner/neighbour
    // for boundary faces: owner
    forAll(faceLabels, faceI)
    {
        // Bug fix: filter inactive faces.  HJ, 21/Mar/2011
        const label faceNo = faceLabels[faceI];

        if (faceNo < nInternalFaces)
        {
            Type t = 0.5*(tf[faceOwner[faceNo]] + tf[faceNeigh[faceNo]]);

            for (direction d = 0; d < nComp; d++)
            {
                vec[d] = component(t, d);
            }
        }
        else if (faceNo < nFaces)
        {
            const Type& t = tf[faceOwner[faceNo]];
            for (direction d = 0; d < nComp; d++)
            {
                vec[d] = component(t, d);
            }
        }
        else
        {
            for (direction d = 0; d < nComp; d++)
            {
                vec[d] = 0;
            }
        }
        vtkFoamTupleRemap<Type>(vec);

        cellData->InsertTuple(faceI, vec);
    }


    vtkPolyData::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    )   ->GetCellData()
        ->AddArray(cellData);

    cellData->Delete();
}


template<class Type>
void Foam::vtkPV4Foam::convertFaceField
(
    const GeometricField<Type, fvPatchField, volMesh>& tf,
    vtkMultiBlockDataSet* output,
    const partInfo& selector,
    const label datasetNo,
    const fvMesh& mesh,
    const faceSet& fSet
)
{
    const label nComp = pTraits<Type>::nComponents;
    const label nInternalFaces = mesh.nInternalFaces();
    const label nFaces = mesh.nFaces();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeigh = mesh.faceNeighbour();

    vtkFloatArray *cellData = vtkFloatArray::New();
    cellData->SetNumberOfTuples( fSet.size() );
    cellData->SetNumberOfComponents( nComp );
    cellData->Allocate( nComp*fSet.size() );
    cellData->SetName( tf.name().c_str() );

    if (debug)
    {
        Info<< "convert convertFaceField: "
            << tf.name()
            << " size = " << tf.size()
            << " nComp=" << nComp
            << " nTuples = " << fSet.size() <<  endl;
    }

    float vec[nComp];

    // for interior faces: average owner/neighbour
    // for boundary faces: owner
    label faceI = 0;
    forAllConstIter(faceSet, fSet, iter)
    {
        const label faceNo = iter.key();

        // Bug fix: filter inactive faces.  HJ, 21/Mar/2011
        if (faceNo < nInternalFaces)
        {
            Type t = 0.5*(tf[faceOwner[faceNo]] + tf[faceNeigh[faceNo]]);

            for (direction d = 0; d < nComp; d++)
            {
                vec[d] = component(t, d);
            }
        }
        else if (faceNo < nFaces)
        {
            const Type& t = tf[faceOwner[faceNo]];
            for (direction d = 0; d < nComp; d++)
            {
                vec[d] = component(t, d);
            }
        }
        else
        {
            for (direction d = 0; d < nComp; d++)
            {
                vec[d] = 0;
            }
        }
        vtkFoamTupleRemap<Type>(vec);

        cellData->InsertTuple(faceI, vec);
        ++faceI;
    }


    vtkPolyData::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    )   ->GetCellData()
        ->AddArray(cellData);

    cellData->Delete();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
