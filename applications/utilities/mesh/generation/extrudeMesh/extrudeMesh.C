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
    Extrude mesh from existing patch or from patch read from file.
    Note: Merges close points so be careful.

    Type of extrusion prescribed by run-time selectable model.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "dimensionedTypes.H"
#include "IFstream.H"
#include "faceMesh.H"
#include "mapPolyMesh.H"
#include "directTopoChange.H"
#include "polyTopoChange.H"
#include "polyTopoChanger.H"
#include "edgeCollapser.H"
#include "mathematicalConstants.H"
#include "globalMeshData.H"
#include "perfectInterface.H"

#include "extrudedMesh.H"
#include "extrudeModel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "setRoots.H"
#   include "createTimeExtruded.H"


    if (args.options().found("sourceRoot") == args.options().found("surface"))
    {
        FatalErrorIn(args.executable())
            << "Need to specify either -sourceRoot/Case/Patch or -surface"
            << " option to specify the source of the patch to extrude"
            << exit(FatalError);
    }


    autoPtr<extrudedMesh> meshPtr(NULL);

    autoPtr<extrudeModel> model
    (
        extrudeModel::New
        (
            IOdictionary
            (
                IOobject
                (
                    "extrudeProperties",
                    runTimeExtruded.constant(),
                    runTimeExtruded,
                    IOobject::MUST_READ
                )
            )
        )
    );

    if (args.options().found("sourceRoot"))
    {
        fileName rootDirSource(args.options()["sourceRoot"]);
        fileName caseDirSource(args.options()["sourceCase"]);
        fileName patchName(args.options()["sourcePatch"]);

        Info<< "Extruding patch " << patchName
            << " on mesh " << rootDirSource << ' ' << caseDirSource << nl
            << endl;

        Time runTime
        (
            Time::controlDictName,
            rootDirSource,
            caseDirSource
        );
#       include "createPolyMesh.H"

        label patchID = mesh.boundaryMesh().findPatchID(patchName);

        if (patchID == -1)
        {
            FatalErrorIn(args.executable())
                << "Cannot find patch " << patchName
                << " in the source mesh.\n"
                << "Valid patch names are " << mesh.boundaryMesh().names()
                << exit(FatalError);
        }

        const polyPatch& pp = mesh.boundaryMesh()[patchID];

        {
            fileName surfName(patchName + ".sMesh");

            Info<< "Writing patch as surfaceMesh to " << surfName << nl << endl;

            faceMesh fMesh(pp.localFaces(), pp.localPoints());

            OFstream os(surfName);
            os << fMesh << nl;
        }

        meshPtr.reset
        (
            new extrudedMesh
            (
                IOobject
                (
                    extrudedMesh::defaultRegion,
                    runTimeExtruded.constant(),
                    runTimeExtruded
                ),
                pp,
                model()
            )
        );
    }
    else
    {
        // Read from surface
        fileName surfName(args.options()["surface"]);

        Info<< "Extruding surfaceMesh read from file " << surfName << nl
            << endl;

        IFstream is(surfName);

        faceMesh fMesh(is);

        Info<< "Read patch from file " << surfName << ':' << nl
            << "    points : " << fMesh.points().size() << nl
            << "    faces  : " << fMesh.size() << nl
            << endl;

        meshPtr.reset
        (
            new extrudedMesh
            (
                IOobject
                (
                    extrudedMesh::defaultRegion,
                    runTimeExtruded.constant(),
                    runTimeExtruded
                ),
                fMesh,
                model()
            )
        );        
    }
    extrudedMesh& mesh = meshPtr();


    const boundBox& bb = mesh.globalData().bb();
    const vector span(bb.max() - bb.min());
    const scalar minDim = min(span[0], min(span[1], span[2]));
    const scalar mergeDim = 1E-4*minDim;

    Pout<< "Mesh bounding box:" << bb << nl
        << "        with span:" << span << nl
        << "Merge distance   :" << mergeDim << nl
        << endl;

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const label origPatchID = patches.findPatchID("originalPatch");
    const label otherPatchID = patches.findPatchID("otherSide");

    if (origPatchID == -1 || otherPatchID == -1)
    {
        FatalErrorIn(args.executable())
            << "Cannot find patch originalPatch or otherSide." << nl
            << "Valid patches are " << patches.names() << exit(FatalError);
    }

    // Collapse edges
    // ~~~~~~~~~~~~~~

    {
        Pout<< "Collapsing edges < " << mergeDim << " ..." << nl << endl;

        // Edge collapsing engine
        edgeCollapser collapser(mesh);

        const edgeList& edges = mesh.edges();
        const pointField& points = mesh.points();

        forAll(edges, edgeI)
        {
            const edge& e = edges[edgeI];

            scalar d = e.mag(points);

            if (d < mergeDim)
            {
                Pout<< "Merging edge " << e << " since length " << d
                    << " << " << mergeDim << nl;

                // Collapse edge to e[0]
                collapser.collapseEdge(edgeI, e[0]);
            }
        }

        // Topo change container
        directTopoChange meshMod(mesh);
        // Put all modifications into meshMod
        bool anyChange = collapser.setRefinement(meshMod);

        if (anyChange)
        {
            // Construct new mesh from directTopoChange.
            autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

            // Update fields
            mesh.updateMesh(map);

            // Move mesh (if inflation used)
            if (map().hasMotionPoints())
            {
                mesh.movePoints(map().preMotionPoints());
            }
        }
    }


    // Merging front and back patch faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (args.options().found("mergeFaces"))
    {
        Pout<< "Assuming full 360 degree axisymmetric case;"
            << " stitching faces on patches " 
            << patches[origPatchID].name() << " and "
            << patches[otherPatchID].name() << " together ..." << nl << endl;

        polyTopoChanger stitcher(mesh);
        stitcher.setSize(1);

        // Make list of masterPatch faces
        labelList isf(patches[origPatchID].size());

        forAll (isf, i)
        {
            isf[i] = patches[origPatchID].start() + i;
        }

        const word cutZoneName("originalCutFaceZone");

        List<faceZone*> fz
        (
            1,
            new faceZone
            (
                cutZoneName,
                isf,
                boolList(isf.size(), false),
                0,
                mesh.faceZones()
            )
        );

        mesh.addZones(List<pointZone*>(0), fz, List<cellZone*>(0));

        // Add the perfect interface mesh modifier
        stitcher.set
        (
            0,
            new perfectInterface
            (
                "couple",
                0,
                stitcher,
                cutZoneName,
                patches[origPatchID].name(),
                patches[otherPatchID].name()
            )
        );

        // Execute all polyMeshModifiers
        autoPtr<mapPolyMesh> morphMap = stitcher.changeMesh();

        mesh.movePoints(morphMap->preMotionPoints());
    }

    if (!mesh.write())
    {
        FatalErrorIn(args.executable()) << "Failed writing mesh"
            << exit(FatalError);
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
