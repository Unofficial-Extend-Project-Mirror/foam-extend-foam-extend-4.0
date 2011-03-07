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
    reconstructPar

Description
    Decompose point, face and cell sets after the case has been decomposed

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "processorMeshes.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
#   include "addRegionOption.H"

#   include "setRootCase.H"
#   include "createTime.H"

    Info<< "Time = " << runTime.timeName() << endl;

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

#   include "createNamedMesh.H"

    // Set all times on processor meshes equal to decomposed mesh
    forAll (databases, procI)
    {
        databases[procI].setTime(runTime.timeName(), runTime.timeIndex());
    }

    // Read all meshes and addressing to reconstructed mesh
    processorMeshes procMeshes(databases, regionName);

    // Find all sets on complete mesh


    // Search for list of objects for the time of the mesh
    IOobjectList objects
    (
        mesh,
        mesh.facesInstance(),
        polyMesh::meshSubDir/"sets"
    );

    Info<< "Searched : " << mesh.facesInstance()/polyMesh::meshSubDir/"sets"
        << nl
        << "Found    : " << objects.names() << nl
        << endl;


    IOobjectList pointObjects(objects.lookupClass(pointSet::typeName));

    for
    (
        IOobjectList::const_iterator iter = pointObjects.begin();
        iter != pointObjects.end();
        ++iter
    )
    {
        // Set not in memory. Load it.
        pointSet set(*iter());

        // Go through all processors
        forAll (procMeshes.meshes(), procI)
        {
            const labelList& addr = procMeshes.pointProcAddressing()[procI];

            labelHashSet procSet;

            forAll (addr, pointI)
            {
                if (set.found(addr[pointI]))
                {
                    procSet.insert(pointI);
                }
            }

            if (!procSet.empty())
            {
                // Set created, write it
                Info<< "Writing point set " << set.name()
                    << " on processor " << procI << endl;
                pointSet cs
                (
                    procMeshes.meshes()[procI],
                    set.name(),
                    procSet,
                    IOobject::NO_WRITE
                );
                cs.write();
            }
        }
    }

    IOobjectList faceObjects(objects.lookupClass(faceSet::typeName));

    for
    (
        IOobjectList::const_iterator iter = faceObjects.begin();
        iter != faceObjects.end();
        ++iter
    )
    {
        // Set not in memory. Load it.
        faceSet set(*iter());

        // Go through all processors
        forAll (procMeshes.meshes(), procI)
        {
            const labelList& addr = procMeshes.faceProcAddressing()[procI];

            labelHashSet procSet;

            forAll (addr, faceI)
            {
                // Note faceProcAddressing peculiarity:
                // change of sign and offset.  HJ, 7/Mar/2011
                if (set.found(mag(addr[faceI]) - 1))
                {
                    procSet.insert(faceI);
                }
            }

            if (!procSet.empty())
            {
                // Set created, write it
                Info<< "Writing face set " << set.name()
                    << " on processor " << procI << endl;
                faceSet cs
                (
                    procMeshes.meshes()[procI],
                    set.name(),
                    procSet,
                    IOobject::NO_WRITE
                );
                cs.write();
            }
        }
    }

    IOobjectList cellObjects(objects.lookupClass(cellSet::typeName));

    for
    (
        IOobjectList::const_iterator iter = cellObjects.begin();
        iter != cellObjects.end();
        ++iter
    )
    {
        // Set not in memory. Load it.
        cellSet set(*iter());

        // Go through all processors
        forAll (procMeshes.meshes(), procI)
        {
            const labelList& addr = procMeshes.cellProcAddressing()[procI];

            labelHashSet procSet;

            forAll (addr, cellI)
            {
                if (set.found(addr[cellI]))
                {
                    procSet.insert(cellI);
                }
            }

            if (!procSet.empty())
            {
                // Set created, write it
                Info<< "Writing cell set " << set.name()
                    << " on processor " << procI << endl;
                cellSet cs
                (
                    procMeshes.meshes()[procI],
                    set.name(),
                    procSet,
                    IOobject::NO_WRITE
                );
                cs.write();
            }
        }
    }


    Info<< "End.\n" << endl;

    return 0;
}


// ************************************************************************* //
