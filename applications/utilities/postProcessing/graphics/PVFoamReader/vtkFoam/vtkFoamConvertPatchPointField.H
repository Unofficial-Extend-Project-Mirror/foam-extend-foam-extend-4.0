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
    Foam::vtkFoam

\*---------------------------------------------------------------------------*/

#ifndef vtkFoamConvertPatchPointField_H
#define vtkFoamConvertPatchPointField_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkFoam::convertPatchPointField
(
    const word& name,
    const Field<Type>& pptf,
    const label regioni
)
{
    vtkUnstructuredGrid *vtkMesh =
        vtkUnstructuredGrid::SafeDownCast(reader_->GetOutput(regioni));

    vtkFloatArray *pointTypes = vtkFloatArray::New();
    pointTypes->SetNumberOfTuples(pptf.size());
    pointTypes->SetNumberOfComponents(Type::nComponents);
    pointTypes->Allocate(Type::nComponents*pptf.size());
    pointTypes->SetName(name.c_str());

    float vec[Type::nComponents];

    forAll(pptf, i)
    {
        for (direction d=0; d<Type::nComponents; d++)
        {
            vec[d] = pptf[i][d];
        }

        pointTypes->InsertTuple(i, vec);
    }

    vtkMesh->GetPointData()->AddArray(pointTypes);
    pointTypes->Delete();
}


template<>
void Foam::vtkFoam::convertPatchPointField
(
    const word& name,
    const Field<scalar>& ppsf,
    const label regioni
)
{
    vtkUnstructuredGrid *vtkMesh =
        vtkUnstructuredGrid::SafeDownCast(reader_->GetOutput(regioni));

    vtkFloatArray *pointScalars = vtkFloatArray::New();
    pointScalars->SetNumberOfTuples(ppsf.size());
    pointScalars->SetNumberOfComponents(1);
    pointScalars->Allocate(ppsf.size());
    pointScalars->SetName(name.c_str());

    for (int i=0; i<ppsf.size(); i++)
    {
        pointScalars->InsertComponent(i, 0, ppsf[i]);
    }

    vtkMesh->GetPointData()->AddArray(pointScalars);
    if (!vtkMesh->GetPointData()->GetScalars())
    {
        vtkMesh->GetPointData()->SetScalars(pointScalars);
    }

    pointScalars->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
