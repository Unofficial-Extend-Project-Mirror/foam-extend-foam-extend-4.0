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
#include "IOobjectList.H"
#include "vtkPV3FoamReader.h"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#include "vtkPV3FoamConvertVolFields.H"
#include "vtkPV3FoamConvertPointFields.H"
#include "vtkPV3FoamConvertLagrangianFields.H"

void Foam::vtkPV3Foam::updateFoamMesh()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateFoamMesh" << endl;
        printMemory();
    }

    if (!reader_->GetCacheMesh())
    {
        delete meshPtr_;
        meshPtr_ = NULL;
    }

    // Check to see if the FOAM mesh has been created
    if (!meshPtr_)
    {
        if (debug)
        {
            Info<< "Creating Foam mesh" << endl;
        }
        meshPtr_ = new fvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                dbPtr_().timeName(),
                dbPtr_()
            )
        );

        meshChanged_ = true;
    }
    else
    {
        if (debug)
        {
            Info<< "Using existing Foam mesh" << endl;
        }
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::updateFoamMesh" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::updateVolFields
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateVolFields" << endl;
        printMemory();
    }

    const fvMesh& mesh = *meshPtr_;

    // Construct interpolation on the raw mesh
    pointMesh pMesh(mesh);

    // Search for list of objects for this time
    IOobjectList objects(mesh, dbPtr_().timeName());

    vtkDataArraySelection* arraySelection = reader_->GetVolFieldSelection();

    // Convert volume fields
    if (debug)
    {
        Info<< "converting Foam volume fields" << endl;
    }

    volPointInterpolation pInterp(mesh, pMesh);

    PtrList<PrimitivePatchInterpolation<primitivePatch> >
        ppInterpList(mesh.boundaryMesh().size());

    forAll(ppInterpList, i)
    {
        ppInterpList.set
        (
            i,
            new PrimitivePatchInterpolation<primitivePatch>
            (
                mesh.boundaryMesh()[i]
            )
        );
    }
/*
    convertVolFields<Foam::label>
    (
        mesh, pInterp, objects, arraySelection, output
    );
*/
    convertVolFields<Foam::scalar>
    (
        mesh, pInterp, ppInterpList, objects, arraySelection, output
    );
    convertVolFields<Foam::vector>
    (
        mesh, pInterp, ppInterpList, objects, arraySelection, output
    );
    convertVolFields<Foam::sphericalTensor>
    (
        mesh, pInterp, ppInterpList, objects, arraySelection, output
    );
    convertVolFields<Foam::symmTensor>
    (
        mesh, pInterp, ppInterpList, objects, arraySelection, output
    );
    convertVolFields<Foam::tensor>
    (
        mesh, pInterp, ppInterpList, objects, arraySelection, output
    );

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::updateVolFields" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::updatePointFields
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updatePointFields" << endl;
        printMemory();
    }

    const fvMesh& mesh = *meshPtr_;

    // Search for list of objects for this time
    IOobjectList objects(mesh, dbPtr_().timeName());

    vtkDataArraySelection* arraySelection = reader_->GetPointFieldSelection();

/*
    convertPointFields<Foam::label>
    (
        mesh, objects, arraySelection, output
    );
*/
    convertPointFields<Foam::scalar>
    (
        mesh, objects, arraySelection, output
    );
    convertPointFields<Foam::vector>
    (
        mesh, objects, arraySelection, output
    );
    convertPointFields<Foam::sphericalTensor>
    (
        mesh, objects, arraySelection, output
    );
    convertPointFields<Foam::symmTensor>
    (
        mesh, objects, arraySelection, output
    );
    convertPointFields<Foam::tensor>
    (
        mesh, objects, arraySelection, output
    );

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::updatePointFields" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::updateLagrangianFields
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateLagrangianFields" << endl;
        printMemory();
    }

    const fvMesh& mesh = *meshPtr_;

    // Search for list of objects for this time
    //- TODO - currently hard-coded to ONE cloud
    IOobjectList objects
    (
        mesh,
        dbPtr_().timeName(),
        "lagrangian"/cloudName_
    );

    vtkDataArraySelection* arraySelection =
        reader_->GetLagrangianFieldSelection();

    // Convert Lagrangian fields
    if (debug)
    {
        Info<< "converting Foam Lagrangian fields" << endl;
    }

    convertLagrangianFields<Foam::label>
    (
        mesh, objects, arraySelection, output
    );

    convertLagrangianFields<Foam::scalar>
    (
        mesh, objects, arraySelection, output
    );
    convertLagrangianFields<Foam::vector>
    (
        mesh, objects, arraySelection, output
    );
    convertLagrangianFields<Foam::sphericalTensor>
    (
        mesh, objects, arraySelection, output
    );
    convertLagrangianFields<Foam::symmTensor>
    (
        mesh, objects, arraySelection, output
    );
    convertLagrangianFields<Foam::tensor>
    (
        mesh, objects, arraySelection, output
    );

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::updateLagrangianFields" << endl;
        printMemory();
    }
}


// ************************************************************************* //
