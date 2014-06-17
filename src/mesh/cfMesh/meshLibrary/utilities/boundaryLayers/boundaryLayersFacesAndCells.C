/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "boundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "demandDrivenData.H"
#include "VRWGraphList.H"

#include <map>

//#define DEBUGLayer

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayers::createNewFacesAndCells(const boolList& treatPatches)
{
    Info << "Starting creating layer cells" << endl;

    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();
    const labelList& faceOwners = mse.faceOwners();

    //- this is used for parallel runs
    const Map<label>* otherProcPatchPtr(NULL);

    if( Pstream::parRun() )
    {
        createNewFacesParallel(treatPatches);

        otherProcPatchPtr = &mse.otherEdgeFacePatch();
    }

    //- create lists for new boundary faces
    VRWGraph newBoundaryFaces;
    labelLongList newBoundaryOwners;
    labelLongList newBoundaryPatches;

    //- create storage for new cells
    VRWGraphList cellsToAdd;

    //- create layer cells and store boundary faces
    const label nOldCells = mesh_.cells().size();
    forAll(bFaces, bfI)
        if( treatPatches[boundaryFacePatches[bfI]] )
        {
            const face& f = bFaces[bfI];

            faceList cellFaces(f.size() + 2);

            label fI(0);

            //- store boundary face
            cellFaces[fI++] = f.reverseFace();

            //- create parallel face
            face newF(f.size());
            forAll(f, pI)
                newF[pI] = newLabelForVertex_[f[pI]];
            cellFaces[fI++] = newF;

            newBoundaryFaces.appendList(newF);
            newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
            newBoundaryPatches.append(boundaryFacePatches[bfI]);

            //- create quad faces
            newF.setSize(4);
            forAll(f, pI)
            {
                newF[0] = f[pI];
                newF[1] = f.nextLabel(pI);
                newF[2] = newLabelForVertex_[f.nextLabel(pI)];
                newF[3] = newLabelForVertex_[f[pI]];

                cellFaces[fI++] = newF;

                //- check if the face is at the boundary
                //- of the treated partitions
                const label edgeI = faceEdges(bfI, pI);
                if( edgeFaces.sizeOfRow(edgeI) == 2 )
                {
                    label neiFace = edgeFaces(edgeI, 0);
                    if( neiFace == bfI )
                        neiFace = edgeFaces(edgeI, 1);

                    if( !treatPatches[boundaryFacePatches[neiFace]] )
                    {
                        newBoundaryFaces.appendList(newF);
                        newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
                        newBoundaryPatches.append(boundaryFacePatches[neiFace]);
                    }
                }
                else if( edgeFaces.sizeOfRow(edgeI) == 1 )
                {
                    const Map<label>& otherProcPatch = *otherProcPatchPtr;
                    if( !treatPatches[otherProcPatch[edgeI]] )
                    {
                        newBoundaryFaces.appendList(newF);
                        newBoundaryOwners.append(cellsToAdd.size() + nOldCells);
                        newBoundaryPatches.append(otherProcPatch[edgeI]);
                    }
                }
            }

            cellsToAdd.appendGraph(cellFaces);
        }
        else
        {
            # ifdef DEBUGLayer
            Info << "Storing original boundary face "
                << bfI << " into patch " << boundaryFacePatches[bfI] << endl;
            # endif

            newBoundaryFaces.appendList(bFaces[bfI]);
            newBoundaryOwners.append(faceOwners[bfI]);
            newBoundaryPatches.append(boundaryFacePatches[bfI]);
        }

    //- create mesh modifier
    polyMeshGenModifier meshModifier(mesh_);

    meshModifier.addCells(cellsToAdd);
    cellsToAdd.clear();
    meshModifier.reorderBoundaryFaces();
    meshModifier.replaceBoundary
    (
        patchNames_,
        newBoundaryFaces,
        newBoundaryOwners,
        newBoundaryPatches
    );

    //- delete meshSurfaceEngine
    this->clearOut();

    # ifdef DEBUGLayer
    mesh_.addressingData().checkMesh(true);
    # endif

    Info << "Finished creating layer cells" << endl;
}

void boundaryLayers::createNewFacesParallel
(
    const boolList& treatPatches
)
{
    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelList& boundaryFacePatches = mse.boundaryFacePatches();
    const labelList& globalEdgeLabel = mse.globalBoundaryEdgeLabel();
    const Map<label>& globalToLocal = mse.globalToLocalBndEdgeAddressing();

    const Map<label>& otherProcPatch = mse.otherEdgeFacePatch();
    const Map<label>& otherFaceProc = mse.otherEdgeFaceAtProc();

    //- the next stage is the generation of processor faces
    //- this step can be done without any communication, but only if the faces
    //- are added in the same order on both processors
    //- this will be achieved by sorting edges according to their global labes
    //- another difficulty here is that new processor patches may occur
    //- during this procedure
    Map<label> otherProcToProcPatch;
    forAll(mesh_.procBoundaries(), patchI)
    {
        const processorBoundaryPatch& wp = mesh_.procBoundaries()[patchI];
        otherProcToProcPatch.insert(wp.neiProcNo(), patchI);
    }

    label nTreatedEdges(0);
    boolList treatEdge(edgeFaces.size(), false);
    for
    (
        Map<label>::const_iterator iter=globalToLocal.begin();
        iter!=globalToLocal.end();
        ++iter
    )
    {
        const label beI = iter();

        if( edgeFaces.sizeOfRow(beI) != 1 )
            continue;

        if(
            treatPatches[boundaryFacePatches[edgeFaces(beI, 0)]] &&
            treatPatches[otherProcPatch[beI]]
        )
        {
            ++nTreatedEdges;
            treatEdge[beI] = true;
        }
    }

    //- create a list of treated edges and sort the list
    labelList treatedEdgeLabels(nTreatedEdges);
    nTreatedEdges = 0;
    forAll(treatEdge, beI)
        if( treatEdge[beI] )
        {
            treatedEdgeLabels[nTreatedEdges++] = globalEdgeLabel[beI];
        }
    treatedEdgeLabels.setSize(nTreatedEdges);

    sort(treatedEdgeLabels);

    //- create additional processor patches if needed
    forAll(treatedEdgeLabels, eI)
    {
        const label beI = globalToLocal[treatedEdgeLabels[eI]];

        if( !otherProcToProcPatch.found(otherFaceProc[beI]) )
        {
            otherProcToProcPatch.insert
            (
                otherFaceProc[beI],
                polyMeshGenModifier(mesh_).addProcessorPatch
                (
                    otherFaceProc[beI]
                )
            );
        }
    }

    //- create new processor faces
    VRWGraph newProcFaces;
    labelLongList faceProcPatch;
    FixedList<label, 4> newF;
    forAll(treatedEdgeLabels, geI)
    {
        const label beI = globalToLocal[treatedEdgeLabels[geI]];

        if( edgeFaces.sizeOfRow(beI) == 0 )
            continue;

        const label bfI = edgeFaces(beI, 0);
        const label pos = faceEdges.containsAtPosition(bfI, beI);
        const edge e = bFaces[bfI].faceEdge(pos);

        if( otherFaceProc[beI] > Pstream::myProcNo() )
        {
            newF[0] = e.start();
            newF[1] = e.end();
            if( patchKey_.size() != 0 )
            {
                newF[2] =
                    findNewNodeLabel(e.end(), patchKey_[otherProcPatch[beI]]);
                newF[3] =
                    findNewNodeLabel(e.start(), patchKey_[otherProcPatch[beI]]);
            }
            else
            {
                newF[2] = newLabelForVertex_[e.end()];
                newF[3] = newLabelForVertex_[e.start()];
            }
        }
        else
        {
            newF[0] = e.end();
            if( patchKey_.size() != 0 )
            {
                newF[1] =
                    findNewNodeLabel
                    (
                        e.end(),
                        patchKey_[boundaryFacePatches[bfI]]
                    );
                newF[2] =
                    findNewNodeLabel
                    (
                        e.start(),
                        patchKey_[boundaryFacePatches[bfI]]
                    );
            }
            else
            {
                newF[1] = newLabelForVertex_[e.end()];
                newF[2] = newLabelForVertex_[e.start()];
            }
            newF[3] = e.start();
        }

        newProcFaces.appendList(newF);
        faceProcPatch.append(otherProcToProcPatch[otherFaceProc[beI]]);
    }

    //- add faces into the mesh
    polyMeshGenModifier(mesh_).addProcessorFaces(newProcFaces, faceProcPatch);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
