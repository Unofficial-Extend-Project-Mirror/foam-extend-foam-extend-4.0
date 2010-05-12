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
    Makes internal faces into boundary faces. Does not duplicate points. Use
    mergeOrSplitBaffles if you want this.

    Note: if any coupled patch face is selected for baffling automatically
    the opposite member is selected for baffling as well. Note that this
    is the same as repatching. This was added only for convenience so
    you don't have to filter coupled boundary out of your set.

\*---------------------------------------------------------------------------*/

#include "syncTools.H"
#include "argList.H"
#include "Time.H"
#include "faceSet.H"
#include "directTopoChange.H"
#include "polyModifyFace.H"
#include "polyAddFace.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("set");
    argList::validArgs.append("patch");
    argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const faceZoneMesh& faceZones = mesh.faceZones();

    // Faces to baffle
    word setName(args.additionalArgs()[0]);
    Pout<< "Reading faceSet from " << setName << nl << endl;
    faceSet facesToSplit(mesh, setName);
    Pout<< "Read " << facesToSplit.size() << " faces from " << setName
        << nl << endl;

    // Patch to put them into
    word patchName(args.additionalArgs()[1]);
    label wantedPatchI = patches.findPatchID(patchName);

    Pout<< "Using patch " << patchName << " at index " << wantedPatchI << endl;

    if (wantedPatchI == -1)
    {
        FatalErrorIn(args.executable())
            << "Cannot find patch " << patchName << exit(FatalError);
    }

    bool overwrite = args.options().found("overwrite");

    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.

    PtrList<volScalarField> vsFlds;
    ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    ReadFields(mesh, objects, vvFlds);

    PtrList<volSphericalTensorField> vstFlds;
    ReadFields(mesh, objects, vstFlds);

    PtrList<volSymmTensorField> vsymtFlds;
    ReadFields(mesh, objects, vsymtFlds);

    PtrList<volTensorField> vtFlds;
    ReadFields(mesh, objects, vtFlds);

    // Read surface fields.

    PtrList<surfaceScalarField> ssFlds;
    ReadFields(mesh, objects, ssFlds);

    PtrList<surfaceVectorField> svFlds;
    ReadFields(mesh, objects, svFlds);

    PtrList<surfaceSphericalTensorField> sstFlds;
    ReadFields(mesh, objects, sstFlds);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    ReadFields(mesh, objects, ssymtFlds);

    PtrList<surfaceTensorField> stFlds;
    ReadFields(mesh, objects, stFlds);


    // Mesh change container
    directTopoChange meshMod(mesh);


    // Creating baffles:
    // - coupled boundary faces : become the patch specified
    // - non-coupled ,,         : illegal
    // - internal faces         : converted into boundary faces.

    labelList newPatch(mesh.nFaces(), -1);

    forAllConstIter(faceSet, facesToSplit, iter)
    {
        label faceI = iter.key();

        label patchI = patches.whichPatch(faceI);

        if (patchI == -1)
        {
            newPatch[faceI] = wantedPatchI;
        }
        else
        {
            if (patches[patchI].coupled())
            {
                if (patchI != wantedPatchI)
                {
                    newPatch[faceI] = wantedPatchI;
                }
            }
            else
            {
                FatalErrorIn(args.executable())
                    << "Can only create baffles from internal faces"
                    << " or coupled boundary faces." << endl
                    << "Face " << faceI << " is a boundary face on patch "
                    << patches[patchI].name() << exit(FatalError);
            }
        }
    }


    // If one side of a coupled boundary is marked for baffling, make sure to
    // also do the other side.

    syncTools::syncFaceList(mesh, newPatch, maxEqOp<label>(), false);


    label nBaffled = 0;

    forAll(newPatch, faceI)
    {
        if (newPatch[faceI] != -1)
        {
            const face& f = mesh.faces()[faceI];
            label zoneID = faceZones.whichZone(faceI);
            bool zoneFlip = false;
            if (zoneID >= 0)
            {
                const faceZone& fZone = faceZones[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    f,                          // modified face
                    faceI,                      // label of face
                    mesh.faceOwner()[faceI],    // owner
                    -1,                         // neighbour
                    false,                      // face flip
                    newPatch[faceI],            // patch for face
                    false,                      // remove from zone
                    zoneID,                     // zone for face
                    zoneFlip                    // face flip in zone
                )
            );

            if (mesh.isInternalFace(faceI))
            {
                meshMod.setAction
                (
                    polyAddFace
                    (
                        f.reverseFace(),            // modified face
                        mesh.faceNeighbour()[faceI],// owner
                        -1,                         // neighbour
                        -1,                         // masterPointID
                        -1,                         // masterEdgeID
                        faceI,                      // masterFaceID,
                        false,                      // face flip
                        newPatch[faceI],            // patch for face
                        zoneID,                     // zone for face
                        zoneFlip                    // face flip in zone
                    )
                );
            }

            nBaffled++;
        }
    }


    Pout<< "Converted locally " << nBaffled
        << " faces into boundary faces on patch " << patchName << nl << endl;

    if (!overwrite)
    {
        runTime++;
    }

    // Change the mesh. Change points directly (no inflation).
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

    // Update fields
    mesh.updateMesh(map);

    // Move mesh (since morphing might not do this)
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }

    Pout<< "Writing mesh to " << runTime.timeName() << endl;

    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
