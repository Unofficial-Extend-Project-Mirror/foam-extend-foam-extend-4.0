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

Application
    reconstructParMesh

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Description
    Reconstructs a mesh using geometrical matching and catenation.
    Use following topological changes in parallel to create global mesh
    and xxxxProcAddressing files in the processor meshes.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "processorMeshesReconstructor.H"
#include "fvFieldReconstructor.H"
#include "pointFieldReconstructor.H"
#include "tetPointFieldReconstructor.H"
#include "reconstructLagrangian.H"

#include "faCFD.H"
#include "faMesh.H"
#include "processorFaMeshes.H"
#include "faFieldReconstructor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // enable -constant ... if someone really wants it
    // enable -zeroTime to prevent accidentally trashing the initial fields
    timeSelector::addOptions(false, true);
    argList::noParallel();
#   include "addRegionOption.H"
    argList::validOptions.insert("cellDist", "");
    argList::validOptions.insert("fields", "\"(list of fields)\"");
    argList::validOptions.insert("noLagrangian", "");

#   include "setRootCase.H"

    bool writeCellDist = args.optionFound("cellDist");

#   include "createTime.H"

    HashSet<word> selectedFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedFields;
    }

    bool noLagrangian = args.optionFound("noLagrangian");

    // Determine the processor count directly
    label nProcs = 0;
    while (isDir(args.path()/(word("processor") + name(nProcs))))
    {
        ++nProcs;
    }

    if (!nProcs)
    {
        FatalErrorIn(args.executable())
            << "No processor* directories found"
            << exit(FatalError);
    }

    // Create the processor databases
    PtrList<Time> databases(nProcs);

    forAll (databases, procI)
    {
        databases.set
        (
            procI,
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/fileName(word("processor") + name(procI))
            )
        );
    }

    // use the times list from the master processor
    // and select a subset based on the command-line options
    instantList timeDirs = timeSelector::select
    (
        databases[0].times(),
        args
    );

    if (timeDirs.empty())
    {
        FatalErrorIn(args.executable())
            << "No times selected"
            << exit(FatalError);
    }

    Foam::word regionName = polyMesh::defaultRegion;

    if (args.optionReadIfPresent("region", regionName))
    {
        Info<< "Selecting region " << regionName << " for time = "
            << runTime.timeName() << Foam::nl << Foam::endl;
    }

    // Set all times on processor meshes equal to reconstructed mesh
    forAll (databases, procI)
    {
        Info<< "Reading database for processor " << procI << endl;

        databases[procI].setTime(runTime.timeName(), runTime.timeIndex());
    }

    // Read all meshes and addressing to reconstructed mesh
    processorMeshesReconstructor procMeshes(databases, regionName);

    autoPtr<fvMesh> meshPtr = procMeshes.reconstructMesh(runTime);


    // Mesh write will be controlled by hand
    meshPtr->write();
    meshPtr->setMotionWriteOpt(IOobject::NO_WRITE);
    meshPtr->setTopoWriteOpt(IOobject::NO_WRITE);

    // Get region prefix for lagrangian
    fileName regionPrefix = "";
    if (regionName != fvMesh::defaultRegion)
    {
        regionPrefix = regionName;
    }


    // Loop over all times
    forAll (timeDirs, timeI)
    {
        // Set time for global database
        runTime.setTime(timeDirs[timeI], timeI);

        Info << "Time = " << runTime.timeName() << endl << endl;

        // Set time for all databases
        forAll (databases, procI)
        {
            databases[procI].setTime(timeDirs[timeI], timeI);
        }

        polyMesh::readUpdateState procStat = procMeshes.readUpdate();

        if (procStat == polyMesh::UNCHANGED)
        {
            Info<< "Mesh unchanged" << endl;

            meshPtr->setMotionWriteOpt(IOobject::NO_WRITE);
            meshPtr->setTopoWriteOpt(IOobject::NO_WRITE);
        }
        else if (procStat == polyMesh::POINTS_MOVED)
        {
            Info<< "Mesh motion detected.  Reconstruct motion points"
                << endl;

            // Reconstruct the points for moving mesh cases and write them out
            procMeshes.reconstructPoints(meshPtr());

            // Set write options
            meshPtr->setMotionWriteOpt(IOobject::AUTO_WRITE);
            meshPtr->setTopoWriteOpt(IOobject::NO_WRITE);

            // Global mesh write
            meshPtr->write();
        }
        else if
        (
            procStat == polyMesh::TOPO_CHANGE
         || procStat == polyMesh::TOPO_PATCH_CHANGE
        )
        {
            Info<< "Topological change detected.  Reconstructing mesh"
                << endl;

            // Reconstruct mesh
            meshPtr = procMeshes.reconstructMesh(runTime);

            // Set write options
            meshPtr->setMotionWriteOpt(IOobject::AUTO_WRITE);
            meshPtr->setTopoWriteOpt(IOobject::AUTO_WRITE);
            procMeshes.writeAddressing();

            // Global mesh write
            meshPtr->write();

            // Write out mapping in processor directories
            forAll (databases, procI)
            {
                databases[procI].write();
            }

            if (writeCellDist)
            {
                // Write as volScalarField for postprocessing.
                volScalarField cellDist
                (
                    IOobject
                    (
                        "cellDist",
                        runTime.timeName(),
                        meshPtr(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    meshPtr(),
                    dimensionedScalar("cellDist", dimless, 0),
                    zeroGradientFvPatchScalarField::typeName
                );
                scalarField& cellDistIn = cellDist.internalField();

                label cellI = 0;

                forAll (procMeshes.meshes(), procI)
                {
                    for
                    (
                        label i = 0;
                        i < procMeshes.meshes()[procI].nCells();
                        i++
                    )
                    {
                        cellDistIn[cellI] = procI;
                        cellI++;
                    }
                }

                cellDist.write();
            }
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Unknown readUpdate state"
                << abort(FatalError);
        }

        fvMesh& mesh = meshPtr();

        // Get list of objects from processor0 database
        IOobjectList objects(procMeshes.meshes()[0], databases[0].timeName());


        // If there are any FV fields, reconstruct them

        if
        (
            objects.lookupClass(volScalarField::typeName).size()
         || objects.lookupClass(volVectorField::typeName).size()
         || objects.lookupClass(volSphericalTensorField::typeName).size()
         || objects.lookupClass(volSymmTensorField::typeName).size()
         || objects.lookupClass(volTensorField::typeName).size()
         || objects.lookupClass(surfaceScalarField::typeName).size()
         || objects.lookupClass(surfaceVectorField::typeName).size()
         || objects.lookupClass(surfaceSphericalTensorField::typeName).size()
         || objects.lookupClass(surfaceSymmTensorField::typeName).size()
         || objects.lookupClass(surfaceTensorField::typeName).size()
        )
        {
            Info << "Reconstructing FV fields" << nl << endl;

            fvFieldReconstructor fvReconstructor
            (
                mesh,
                procMeshes.meshes(),
                procMeshes.faceProcAddressing(),
                procMeshes.cellProcAddressing(),
                procMeshes.boundaryProcAddressing()
            );

            fvReconstructor.reconstructFvVolumeFields<scalar>
            (
                objects,
                selectedFields
            );
            fvReconstructor.reconstructFvVolumeFields<vector>
            (
                objects,
                selectedFields
            );
            fvReconstructor.reconstructFvVolumeFields<sphericalTensor>
            (
                objects,
                selectedFields
            );
            fvReconstructor.reconstructFvVolumeFields<symmTensor>
            (
                objects,
                selectedFields
            );
            fvReconstructor.reconstructFvVolumeFields<tensor>
            (
                objects,
                selectedFields
            );

            fvReconstructor.reconstructFvSurfaceFields<scalar>
            (
                objects,
                selectedFields
            );
            fvReconstructor.reconstructFvSurfaceFields<vector>
            (
                objects,
                selectedFields
            );
            fvReconstructor.reconstructFvSurfaceFields<sphericalTensor>
            (
                objects,
                selectedFields
            );
            fvReconstructor.reconstructFvSurfaceFields<symmTensor>
            (
                objects,
                selectedFields
            );
            fvReconstructor.reconstructFvSurfaceFields<tensor>
            (
                objects,
                selectedFields
            );
        }
        else
        {
            Info << "No FV fields" << nl << endl;
        }


        // If there are any point fields, reconstruct them
        if
        (
            objects.lookupClass(pointScalarField::typeName).size()
         || objects.lookupClass(pointVectorField::typeName).size()
         || objects.lookupClass(pointSphericalTensorField::typeName).size()
         || objects.lookupClass(pointSymmTensorField::typeName).size()
         || objects.lookupClass(pointTensorField::typeName).size()
        )
        {
            Info << "Reconstructing point fields" << nl << endl;

            pointMesh pMesh(mesh);
            PtrList<pointMesh> pMeshes(procMeshes.meshes().size());

            forAll (pMeshes, procI)
            {
                pMeshes.set(procI, new pointMesh(procMeshes.meshes()[procI]));
            }

            pointFieldReconstructor pointReconstructor
            (
                pMesh,
                pMeshes,
                procMeshes.pointProcAddressing(),
                procMeshes.boundaryProcAddressing()
            );

            pointReconstructor.reconstructFields<scalar>(objects);
            pointReconstructor.reconstructFields<vector>(objects);
            pointReconstructor.reconstructFields<sphericalTensor>(objects);
            pointReconstructor.reconstructFields<symmTensor>(objects);
            pointReconstructor.reconstructFields<tensor>(objects);
        }
        else
        {
            Info << "No point fields" << nl << endl;
        }

        // If there are any tetFem fields, reconstruct them
        if
        (
            objects.lookupClass(tetPointScalarField::typeName).size()
         || objects.lookupClass(tetPointVectorField::typeName).size()
         || objects.lookupClass(tetPointSphericalTensorField::typeName).size()
         || objects.lookupClass(tetPointSymmTensorField::typeName).size()
         || objects.lookupClass(tetPointTensorField::typeName).size()

         || objects.lookupClass(elementScalarField::typeName).size()
         || objects.lookupClass(elementVectorField::typeName).size()
        )
        {
            Info << "Reconstructing tet point fields" << nl << endl;

            tetPolyMesh tetMesh(mesh);
            PtrList<tetPolyMesh> tetMeshes(procMeshes.meshes().size());

            forAll (tetMeshes, procI)
            {
                tetMeshes.set
                (
                    procI,
                    new tetPolyMesh(procMeshes.meshes()[procI])
                );
            }

            tetPointFieldReconstructor tetPointReconstructor
            (
                tetMesh,
                tetMeshes,
                procMeshes.pointProcAddressing(),
                procMeshes.faceProcAddressing(),
                procMeshes.cellProcAddressing(),
                procMeshes.boundaryProcAddressing()
            );

            // Reconstruct tet point fields
            tetPointReconstructor.reconstructTetPointFields<scalar>(objects);
            tetPointReconstructor.reconstructTetPointFields<vector>(objects);
            tetPointReconstructor.
                reconstructTetPointFields<sphericalTensor>(objects);
            tetPointReconstructor.
                reconstructTetPointFields<symmTensor>(objects);
            tetPointReconstructor.reconstructTetPointFields<tensor>(objects);

            tetPointReconstructor.reconstructElementFields<scalar>(objects);
            tetPointReconstructor.reconstructElementFields<vector>(objects);
        }
        else
        {
            Info << "No tetFem fields" << nl << endl;
        }


        // If there are any clouds, reconstruct them.
        // The problem is that a cloud of size zero will not get written so
        // in pass 1 we determine the cloud names and per cloud name the
        // fields. Note that the fields are stored as IOobjectList from
        // the first processor that has them. They are in pass2 only used
        // for name and type (scalar, vector etc).

        if (!noLagrangian)
        {
            HashTable<IOobjectList> cloudObjects;

            forAll (databases, procI)
            {
                fileNameList cloudDirs
                (
                    readDir
                    (
                        databases[procI].timePath()/regionPrefix/cloud::prefix,
                        fileName::DIRECTORY
                    )
                );

                forAll (cloudDirs, i)
                {
                    // Check if we already have cloud objects for
                    // this cloudname
                    HashTable<IOobjectList>::const_iterator iter =
                        cloudObjects.find(cloudDirs[i]);

                    if (iter == cloudObjects.end())
                    {
                        // Do local scan for valid cloud objects
                        IOobjectList sprayObjs
                        (
                            procMeshes.meshes()[procI],
                            databases[procI].timeName(),
                            cloud::prefix/cloudDirs[i]
                        );

                        IOobject* positionsPtr = sprayObjs.lookup("positions");

                        if (positionsPtr)
                        {
                            cloudObjects.insert(cloudDirs[i], sprayObjs);
                        }
                    }
                }
            }


            if (cloudObjects.size())
            {
                // Pass2: reconstruct the cloud
                forAllConstIter(HashTable<IOobjectList>, cloudObjects, iter)
                {
                    const word cloudName = string::validate<word>(iter.key());

                    // Objects (on arbitrary processor)
                    const IOobjectList& sprayObjs = iter();

                    Info<< "Reconstructing lagrangian fields for cloud "
                        << cloudName << nl << endl;

                    reconstructLagrangianPositions
                    (
                        mesh,
                        cloudName,
                        procMeshes.meshes(),
                        procMeshes.faceProcAddressing(),
                        procMeshes.cellProcAddressing()
                    );
                    reconstructLagrangianFields<label>
                    (
                        cloudName,
                        mesh,
                        procMeshes.meshes(),
                        sprayObjs
                    );
                    reconstructLagrangianFields<scalar>
                    (
                        cloudName,
                        mesh,
                        procMeshes.meshes(),
                        sprayObjs
                    );
                    reconstructLagrangianFields<vector>
                    (
                        cloudName,
                        mesh,
                        procMeshes.meshes(),
                        sprayObjs
                    );
                    reconstructLagrangianFields<sphericalTensor>
                    (
                        cloudName,
                        mesh,
                        procMeshes.meshes(),
                        sprayObjs
                    );
                    reconstructLagrangianFields<symmTensor>
                    (
                        cloudName,
                        mesh,
                        procMeshes.meshes(),
                        sprayObjs
                    );
                    reconstructLagrangianFields<tensor>
                    (
                        cloudName,
                        mesh,
                        procMeshes.meshes(),
                        sprayObjs
                    );
                }
            }
            else
            {
                Info << "No lagrangian fields" << nl << endl;
            }
        }

        // If there are any FA fields, reconstruct them

        if
        (
            objects.lookupClass(areaScalarField::typeName).size()
         || objects.lookupClass(areaVectorField::typeName).size()
         || objects.lookupClass(areaSphericalTensorField::typeName).size()
         || objects.lookupClass(areaSymmTensorField::typeName).size()
         || objects.lookupClass(areaTensorField::typeName).size()
         || objects.lookupClass(edgeScalarField::typeName).size()
        )
        {
            Info << "Reconstructing FA fields" << nl << endl;

            faMesh aMesh(mesh);

            processorFaMeshes procFaMeshes(procMeshes.meshes());

            faFieldReconstructor faReconstructor
            (
                aMesh,
                procFaMeshes.meshes(),
                procFaMeshes.edgeProcAddressing(),
                procFaMeshes.faceProcAddressing(),
                procFaMeshes.boundaryProcAddressing()
            );

            faReconstructor.reconstructFaAreaFields<scalar>(objects);
            faReconstructor.reconstructFaAreaFields<vector>(objects);
            faReconstructor
               .reconstructFaAreaFields<sphericalTensor>(objects);
            faReconstructor.reconstructFaAreaFields<symmTensor>(objects);
            faReconstructor.reconstructFaAreaFields<tensor>(objects);

            faReconstructor.reconstructFaEdgeFields<scalar>(objects);
        }
        else
        {
            Info << "No FA fields" << nl << endl;
        }

        // If there are any "uniform" directories copy them from
        // the master processor

        fileName uniformDir0 = databases[0].timePath()/"uniform";
        if (isDir(uniformDir0))
        {
            cp(uniformDir0, runTime.timePath());
        }
    }

    Info<< "End.\n" << endl;

    return 0;
}


// ************************************************************************* //
