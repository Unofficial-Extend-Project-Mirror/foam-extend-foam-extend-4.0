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
    makeFaMesh

Description
    A mesh generator for finite area mesh.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "argList.H"
#include "OSspecific.H"
#include "faMesh.H"
#include "fvMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class faPatchData
{
public:
    word name_;
    word type_;
    dictionary dict_;
    label ownPolyPatchID_;
    label ngbPolyPatchID_;
    labelList edgeLabels_;
    faPatchData()
    :
        name_(word::null),
        type_(word::null),
        ownPolyPatchID_(-1),
        ngbPolyPatchID_(-1)
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    // Reading faMeshDefinition dictionary
    IOdictionary faMeshDefinition
    (
        IOobject
        (
            "faMeshDefinition",
            runTime.constant(),
            "faMesh",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    wordList polyMeshPatches
    (
        faMeshDefinition.lookup("polyMeshPatches")
    );

    dictionary bndDict = faMeshDefinition.subDict("boundary");

    wordList faPatchNames = bndDict.toc();

    List<faPatchData> faPatches(faPatchNames.size()+1);

    forAll (faPatchNames, patchI)
    {
        dictionary curPatchDict =
            bndDict.subDict(faPatchNames[patchI]);

        faPatches[patchI].name_ = faPatchNames[patchI];

        faPatches[patchI].type_ = 
            word(curPatchDict.lookup("type"));

        faPatches[patchI].ownPolyPatchID_ =
            mesh.boundaryMesh().findPatchID
            (
                word(curPatchDict.lookup("ownerPolyPatch"))
            );

        faPatches[patchI].ngbPolyPatchID_  =
            mesh.boundaryMesh().findPatchID
            (
                word(curPatchDict.lookup("neighbourPolyPatch"))
            );
    }


    // Setting faceLabels list size
    label size = 0;

    labelList patchIDs(polyMeshPatches.size(), -1);

    forAll (polyMeshPatches, patchI)
    {
        patchIDs[patchI] =
            mesh.boundaryMesh().findPatchID(polyMeshPatches[patchI]);

        size += mesh.boundaryMesh()[patchIDs[patchI]].size();
    }

    labelList faceLabels(size, -1);

    sort(patchIDs);


    // Filling of faceLabels list
    label faceI = -1;

    forAll (polyMeshPatches, patchI)
    {
        label start = mesh.boundaryMesh()[patchIDs[patchI]].start();
        
        label size  = mesh.boundaryMesh()[patchIDs[patchI]].size();

        for(label i = 0; i < size; i++)
        {
            faceLabels[++faceI] = start + i;
        }
    }

    // Creating faMesh
    Info << "Create faMesh ... ";

    faMesh areaMesh
    (
        mesh,
        faceLabels
    );
    Info << "Done" << endl;


    // Determination of faPatch ID for each boundary edge.
    // Result is in the bndEdgeFaPatchIDs list
    const indirectPrimitivePatch& patch = areaMesh.patch();

    labelList faceCells(faceLabels.size(), -1);

    forAll (faceCells, faceI)
    {
        label faceID = faceLabels[faceI];

        faceCells[faceI] = mesh.faceOwner()[faceID];
    }

    labelList meshEdges =
        patch.meshEdges
        (
            mesh.edges(),
            mesh.cellEdges(),
            faceCells
        );

    const labelListList& edgeFaces = mesh.edgeFaces();

    const label nTotalEdges = patch.nEdges();
    const label nInternalEdges = patch.nInternalEdges();

    labelList bndEdgeFaPatchIDs(nTotalEdges - nInternalEdges, -1);

    for (label edgeI = nInternalEdges; edgeI < nTotalEdges; edgeI++)
    {
        label curMeshEdge = meshEdges[edgeI];

        labelList curEdgePatchIDs(2, -1);

        label patchI = -1;

        forAll (edgeFaces[curMeshEdge], faceI)
        {
            label curFace = edgeFaces[curMeshEdge][faceI];

            label curPatchID = mesh.boundaryMesh().whichPatch(curFace);

            if (curPatchID != -1)
            {
                curEdgePatchIDs[++patchI] = curPatchID;
            }
        }

        for(label pI = 0; pI < faPatches.size() - 1; pI++)
        {
            if
            (
                (
                    curEdgePatchIDs[0] == faPatches[pI].ownPolyPatchID_
                 && curEdgePatchIDs[1] == faPatches[pI].ngbPolyPatchID_
                )
                ||
                (
                    curEdgePatchIDs[1] == faPatches[pI].ownPolyPatchID_
                 && curEdgePatchIDs[0] == faPatches[pI].ngbPolyPatchID_
                )
            )
            {
                bndEdgeFaPatchIDs[edgeI - nInternalEdges] = pI;
                break;
            }
        }
    }


    // Set edgeLabels for each faPatch
    for(label pI=0; pI<(faPatches.size()-1); pI++)
    {
        SLList<label> tmpList;

        forAll (bndEdgeFaPatchIDs, eI)
        {
            if (bndEdgeFaPatchIDs[eI] == pI)
            {
                tmpList.append(nInternalEdges + eI);
            }
        }

        faPatches[pI].edgeLabels_ = tmpList;
    }

    // Check for undefined edges
    SLList<label> tmpList;

    forAll (bndEdgeFaPatchIDs, eI)
    {
        if (bndEdgeFaPatchIDs[eI] == -1)
        {
            tmpList.append(nInternalEdges + eI);
        }
    }

    if (tmpList.size() > 0)
    {
        label pI = faPatches.size()-1;

        faPatches[pI].name_ = "undefined";
        faPatches[pI].type_ = "patch";
        faPatches[pI].edgeLabels_ = tmpList;
    }

    // Add good patches to faMesh
    SLList<faPatch*> faPatchLst;

    for(label pI = 0; pI < faPatches.size(); pI++)
    {
        faPatches[pI].dict_.add("type", faPatches[pI].type_);
        faPatches[pI].dict_.add("edgeLabels", faPatches[pI].edgeLabels_);
        faPatches[pI].dict_.add
        (
            "ngbPolyPatchIndex",
            faPatches[pI].ngbPolyPatchID_
        );

        if(faPatches[pI].edgeLabels_.size() > 0)
        {
            faPatchLst.append
            (
                faPatch::New
                (
                    faPatches[pI].name_,
                    faPatches[pI].dict_,
                    pI,
                    areaMesh.boundary()
                ).ptr()
            );
        }
    }

    Info << "Add faPatches ... ";
    areaMesh.addFaPatches(List<faPatch*>(faPatchLst));
    Info << "Done" << endl;

    // Writing faMesh
    Info << "Write finite area mesh ... ";
    areaMesh.write();

    Info << "Done" << endl;

    return(0);
}

// ************************************************************************* //
