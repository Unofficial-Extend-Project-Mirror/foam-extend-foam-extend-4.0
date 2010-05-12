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
#include "Cloud.H"
#include "fvMesh.H"
#include "passiveParticle.H"
#include "vtkPV3FoamInsertNextPoint.H"
#include "IOobjectList.H"

// VTK includes
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPV3Foam::addLagrangianMesh
(
    const fvMesh& mesh,
    vtkPolyData *vtkmesh
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::addLagrangianMesh - timePath "
            << mesh.time().timePath()/"lagrangian" << endl;
        printMemory();
    }

    fileNameList cloudDirs
    (
        readDir(mesh.time().timePath()/"lagrangian", fileName::DIRECTORY)
    );

    if (debug && cloudDirs.size())
    {
        Info<< "... check cloudDirs: " << cloudDirs << endl;
    }

    bool foundCloud = false;
    forAll(cloudDirs, i)
    {
        IOobjectList sprayObjs
        (
            mesh,
            mesh.time().timeName(),
            "lagrangian"/cloudDirs[i]
        );

        IOobject* positionsPtr = sprayObjs.lookup("positions");

        if (positionsPtr && !foundCloud)
        {
            foundCloud = true;

            Cloud<passiveParticle> parcels(mesh, cloudDirs[i], false);

            if (debug)
            {
                Info<< "cloud with " << parcels.size() << " parcels" << endl;
            }

            vtkPoints* vtkpoints = vtkPoints::New();
            vtkpoints->Allocate(parcels.size());

            forAllConstIter(Cloud<passiveParticle>, parcels, elmnt)
            {
                vtkPV3FoamInsertNextPoint(vtkpoints, elmnt().position());
            }

            vtkmesh->SetPoints(vtkpoints);
            vtkpoints->Delete();
        }
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::addLagrangianMesh" << endl;
        printMemory();
    }
}


// ************************************************************************* //
