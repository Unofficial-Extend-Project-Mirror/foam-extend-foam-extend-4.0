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
#include "faceSet.H"
#include "vtkPV3FoamInsertNextPoint.H"

// VTK includes
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPV3Foam::addFaceSetMesh
(
    const fvMesh& mesh,
    const faceSet& fSet,
    vtkPolyData* vtkmesh
)
{
    if (debug)
    {
        Info<< "entered add face set internal mesh" << endl;
    }

    // Construct primitivePatch of faces in fSet.

    const faceList& meshFaces = mesh.allFaces();
    faceList patchFaces(fSet.size());
    label faceI = 0;
    forAllConstIter(faceSet, fSet, iter)
    {
        patchFaces[faceI++] = meshFaces[iter.key()];
    }

    // Use all points: support for inactive points and faces.
    // HJ, 28/Mar/2009
    primitiveFacePatch p(patchFaces, mesh.allPoints());


    // The balance of this routine should be identical to addPatchMesh

    // Convert Foam mesh vertices to VTK
    const pointField& points = p.localPoints();

    vtkPoints *vtkpoints = vtkPoints::New();
    vtkpoints->Allocate(points.size());
    forAll(points, i)
    {
        vtkPV3FoamInsertNextPoint(vtkpoints, points[i]);
    }
    vtkmesh->SetPoints(vtkpoints);
    vtkpoints->Delete();

    // Add faces as polygons
    const faceList& faces = p.localFaces();

    vtkCellArray* vtkcells = vtkCellArray::New();
    vtkcells->Allocate(p.size());
    forAll(faces, faceI)
    {
        const face& f = faces[faceI];
        vtkIdType nodeIds[f.size()];

        forAll (f, fp)
        {
            nodeIds[fp] = f[fp];
        }
        vtkcells->InsertNextCell(f.size(), nodeIds);
    }

    vtkmesh->SetPolys(vtkcells);
    vtkcells->Delete();
}

// ************************************************************************* //
