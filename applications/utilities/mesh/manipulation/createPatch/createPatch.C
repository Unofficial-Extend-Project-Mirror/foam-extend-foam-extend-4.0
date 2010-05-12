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
    Utility to create patches out of selected boundary faces. Faces come either
    from existing patches or from a faceSet.

    More specifically it:
    - creates new patches (from selected boundary faces). Synchronise faces
      on coupled patches.
    - synchronises points on coupled boundaries
    - remove patches with 0 faces in them

\*---------------------------------------------------------------------------*/

#include "syncTools.H"
#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "SortableList.H"
#include "OFstream.H"
#include "meshTools.H"
#include "faceSet.H"
#include "IOPtrList.H"
#include "mapPolyMesh.H"
#include "directTopoChange.H"
#include "polyModifyFace.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<dictionary>, 0);
}

// Combine operator to synchronise points. We choose point nearest to origin so
// we can use e.g. great,great,great as null value.
class nearestEqOp
{

public:

    void operator()(vector& x, const vector& y) const
    {
        if (magSqr(y) < magSqr(x))
        {
            x = y;
        }
    }
};



label getPatch(const polyBoundaryMesh& patches, const word& patchName)
{
    label patchI = patches.findPatchID(patchName);

    if (patchI == -1)
    {
        FatalErrorIn("createPatch(const polyBoundaryMesh&, const word&)")
            << "Cannot find source patch " << patchName
            << endl << "Valid patch names are " << patches.names()
            << exit(FatalError);
    }

    return patchI;
}


void changePatchID
(
    const polyMesh& mesh,
    const label faceID,
    const label patchID,
    directTopoChange& meshMod
)
{
    const label zoneID = mesh.faceZones().whichZone(faceID);

    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceID)];
    }

    meshMod.setAction
    (
        polyModifyFace
        (
            mesh.faces()[faceID],               // face
            faceID,                             // face ID
            mesh.faceOwner()[faceID],           // owner
            -1,                                 // neighbour
            false,                              // flip flux
            patchID,                            // patch ID
            false,                              // remove from zone
            zoneID,                             // zone ID
            zoneFlip                            // zone flip
        )
    );
}


// Filter out the empty patches.
void filterPatches(polyMesh& mesh)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Patches to keep
    DynamicList<polyPatch*> allPatches(patches.size());

    label nOldPatches = returnReduce(patches.size(), sumOp<label>());

    // Copy old patches.
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        // Note: reduce possible since non-proc patches guaranteed in same order
        if (!isA<processorPolyPatch>(pp))
        {
            if (returnReduce(pp.size(), sumOp<label>()) > 0)
            {
                allPatches.append
                (
                    pp.clone
                    (
                        patches,
                        allPatches.size(),
                        pp.size(),
                        pp.start()
                    ).ptr()
                );
            }
            else
            {
                Info<< "Removing empty patch " << pp.name() << " at position "
                    << patchI << endl;
            }
        }
    }
    // Copy non-empty processor patches
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            if (pp.size() > 0)
            {
                allPatches.append
                (
                    pp.clone
                    (
                        patches,
                        allPatches.size(),
                        pp.size(),
                        pp.start()
                    ).ptr()
                );
            }
            else
            {
                Info<< "Removing empty processor patch " << pp.name()
                    << " at position " << patchI << endl;
            }
        }
    }

    label nAllPatches = returnReduce(allPatches.size(), sumOp<label>());
    if (nAllPatches != nOldPatches)
    {
        Info<< "Removing patches." << endl;
        allPatches.shrink();
        mesh.removeBoundary();
        mesh.addPatches(allPatches);
    }
    else
    {
        Info<< "No patches removed." << endl;
    }
}


// Dump for all patches the current match
void dumpCyclicMatch(const fileName& prefix, const polyMesh& mesh)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            label halfSize = cycPatch.size()/2;

            // Dump halves
            {
                OFstream str(prefix+cycPatch.name()+"_half0.obj");
                meshTools::writeOBJ
                (
                    str,
                    static_cast<faceList>
                    (
                        SubList<face>
                        (
                            cycPatch,
                            halfSize
                        )
                    ),
                    cycPatch.points()
                );
            }
            {
                OFstream str(prefix+cycPatch.name()+"_half1.obj");
                meshTools::writeOBJ
                (
                    str,
                    static_cast<faceList>
                    (
                        SubList<face>
                        (
                            cycPatch,
                            halfSize,
                            halfSize
                        )
                    ),
                    cycPatch.points()
                );
            }

//            cycPatch.writeOBJ
//            (
//                prefix+cycPatch.name()+"_half0.obj",
//                SubList<face>
//                (
//                    cycPatch,
//                    halfSize
//                ),
//                cycPatch.points()
//            );
//            cycPatch.writeOBJ
//            (
//                prefix+cycPatch.name()+"_half1.obj",
//                SubList<face>
//                (
//                    cycPatch,
//                    halfSize,
//                    halfSize
//                ),
//                cycPatch.points()
//            );

            // Lines between corresponding face centres
            OFstream str(prefix+cycPatch.name()+"_match.obj");
            label vertI = 0;

            for (label faceI = 0; faceI < halfSize; faceI++)
            {
                const point& fc0 = mesh.faceCentres()[cycPatch.start()+faceI];
                meshTools::writeOBJ(str, fc0);
                vertI++;

                label nbrFaceI = halfSize + faceI;
                const point& fc1 = mesh.faceCentres()[cycPatch.start()+nbrFaceI];
                meshTools::writeOBJ(str, fc1);
                vertI++;

                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }
}


// Main program:

int main(int argc, char *argv[])
{
    argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"
#   include "createTime.H"

    const bool overwrite = args.options().found("overwrite");

    Info<< "Reading createPatchDict." << nl << endl;

    IOdictionary dict
    (
        IOobject
        (
            "createPatchDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );


    // Whether to synchronise points
    const Switch pointSync(dict.lookup("pointSync"));


    // Set the matching tolerance so we can read illegal meshes
    scalar tol = readScalar(dict.lookup("matchTolerance"));
    Info<< "Using relative tolerance " << tol
        << " to match up faces and points" << nl << endl;
    // Change tolerancein controlDict instead.  HJ, 22/Oct/2008
//     polyPatch::matchTol_ = tol;


#   include "createPolyMesh.H"

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // If running parallel check same patches everywhere
    patches.checkParallelSync(true);


    dumpCyclicMatch("initial_", mesh);

    // Read patch construct info from dictionary
    PtrList<dictionary> patchSources(dict.lookup("patches"));



    // 1. Add all new patches
    // ~~~~~~~~~~~~~~~~~~~~~~

    if (patchSources.size() > 0)
    {
        // Old and new patches.
        DynamicList<polyPatch*> allPatches(patches.size()+patchSources.size());

        label startFaceI = mesh.nInternalFaces();

        // Copy old patches.
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (!isA<processorPolyPatch>(pp))
            {
                allPatches.append
                (
                    pp.clone
                    (
                        patches,
                        patchI,
                        pp.size(),
                        startFaceI
                    ).ptr()
                );
                startFaceI += pp.size();
            }
        }

        forAll(patchSources, addedI)
        {
            const dictionary& dict = patchSources[addedI];

            word patchName(dict.lookup("name"));

            label destPatchI = patches.findPatchID(patchName);

            word patchType(dict.lookup("type"));

            if (destPatchI == -1)
            {
                destPatchI = allPatches.size();

                Info<< "Adding new patch " << patchName
                    << " of type " << patchType
                    << " as patch " << destPatchI << endl;

                // Add an empty patch.
                allPatches.append
                (
                    polyPatch::New
                    (
                        patchType,
                        patchName,
                        0,              // size
                        startFaceI,     // start
                        destPatchI,
                        patches
                    ).ptr()
                );
            }
        }

        // Copy old patches.
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (isA<processorPolyPatch>(pp))
            {
                allPatches.append
                (
                    pp.clone
                    (
                        patches,
                        patchI,
                        pp.size(),
                        startFaceI
                    ).ptr()
                );
                startFaceI += pp.size();
            }
        }

        allPatches.shrink();
        mesh.removeBoundary();
        mesh.addPatches(allPatches);

        Info<< endl;
    }



    // 2. Repatch faces
    // ~~~~~~~~~~~~~~~~

    directTopoChange meshMod(mesh);


    forAll(patchSources, addedI)
    {
        const dictionary& dict = patchSources[addedI];

        word patchName(dict.lookup("name"));

        label destPatchI = patches.findPatchID(patchName);

        if (destPatchI == -1)
        {
            FatalErrorIn(args.executable()) << "patch " << patchName
                << " not added. Problem." << abort(FatalError);
        }

        word sourceType(dict.lookup("constructFrom"));

        if (sourceType == "patches")
        {
            wordList patchSources(dict.lookup("patches"));

            // Repatch faces of the patches.
            forAll(patchSources, sourceI)
            {
                label patchI = getPatch(patches, patchSources[sourceI]);

                const polyPatch& pp = patches[patchI];

                Info<< "Moving faces from patch " << pp.name()
                    << " to patch " << destPatchI << endl;

                forAll(pp, i)
                {
                    changePatchID
                    (
                        mesh,
                        pp.start() + i,
                        destPatchI,
                        meshMod
                    );
                }
            }
        }
        else if (sourceType == "set")
        {
            word setName(dict.lookup("set"));

            faceSet faces(mesh, setName);

            Info<< "Read " << returnReduce(faces.size(), sumOp<label>())
                << " faces from faceSet " << faces.name() << endl;

            // Sort (since faceSet contains faces in arbitrary order)
            labelList faceLabels(faces.toc());

            SortableList<label> patchFaces(faceLabels);

            forAll(patchFaces, i)
            {
                label faceI = patchFaces[i];

                if (mesh.isInternalFace(faceI))
                {
                    FatalErrorIn(args.executable())
                        << "Face " << faceI << " specified in set "
                        << faces.name()
                        << " is not an external face of the mesh." << endl
                        << "This application can only repatch existing boundary"
                        << " faces." << exit(FatalError);
                }

                changePatchID
                (
                    mesh,
                    faceI,
                    destPatchI,
                    meshMod
                );
            }
        }
        else
        {
            FatalErrorIn(args.executable())
                << "Invalid source type " << sourceType << endl
                << "Valid source types are 'patches' 'set'" << exit(FatalError);
        }
    }
    Info<< endl;


    // Change mesh, use inflation to reforce calculation of transformation
    // tensors.
    Info<< "Doing topology modification to order faces." << nl << endl;
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, true);
    mesh.movePoints(map().preMotionPoints());

    // Synchronise points.
    if (!pointSync)
    {
        Info<< "Not synchronising points." << nl << endl;
    }
    else
    {
        Info<< "Synchronising points." << endl;

        pointField newPoints(mesh.points());
        syncTools::syncPointList
        (
            mesh,
            newPoints,
            nearestEqOp(),              // cop
            point(GREAT, GREAT, GREAT), // nullValue
            true                        // applySeparation
        );

        scalarField diff(mag(newPoints-mesh.points()));
        Info<< "Points changed by average:" << gAverage(diff)
            << " max:" << gMax(diff) << nl << endl;

        mesh.movePoints(newPoints);
    }

    // 3. Remove zeros-sized patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Removing patches with no faces in them." << nl<< endl;
    filterPatches(mesh);



    // Set the precision of the points data to 10
    IOstream::defaultPrecision(10);

    if (!overwrite)
    {
        runTime++;
    }

    // Write resulting mesh
    Info<< "Writing repatched mesh to " << runTime.timeName() << nl << endl;
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
