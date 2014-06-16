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

#include "cartesianMeshExtractor.H"
#include "demandDrivenData.H"
#include "meshOctree.H"
#include "labelledPair.H"

#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "helperFunctionsPar.H"
#include "decomposeFaces.H"
#include "decomposeCells.H"

#include <map>
#include <sstream>

//#define DEBUGMesh

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void cartesianMeshExtractor::createPolyMesh()
{
    Info << "Creating polyMesh from octree" << endl;

    const meshOctree& octree = octreeCheck_.octree();

    //- give labels to cubes which will be used as mesh cells
    const List<direction>& cType = octreeCheck_.boxType();

    labelList& leafCellLabel = *leafCellLabelPtr_;
    label nCells(0);
    forAll(cType, leafI)
    {
        if(
            Pstream::parRun() &&
            (octree.returnLeaf(leafI).procNo() != Pstream::myProcNo())
        )
            continue;

        if( cType[leafI] & meshOctreeAddressing::MESHCELL )
        {
            leafCellLabel[leafI] = nCells++;
        }
    }

    //- access to mesh data
    polyMeshGenModifier meshModifier(mesh_);
    faceListPMG& faces = meshModifier.facesAccess();
    cellListPMG& cells = meshModifier.cellsAccess();

    //- start creating octree mesh
    cells.setSize(nCells);
    List<direction> nFacesInCell(nCells, direction(0));
    label nFaces(0);

    const VRWGraph& octreeFaces = octreeCheck_.octreeFaces();
    const labelLongList& owner = octreeCheck_.octreeFaceOwner();
    const labelLongList& neighbour = octreeCheck_.octreeFaceNeighbour();

    //- map storing box label and a direction for each processor face
    //- The map stores data in the same order on both sides of processor
    //- boundaries. This is a consequence of Morton ordering of
    //- leaf boxes in the octree.
    std::map<label, labelLongList> procFaces;

    forAll(octreeFaces, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        const label ownLabel = leafCellLabel[own];
        label neiLabel(-1);
        if( nei != -1 )
            neiLabel = leafCellLabel[nei];

        if( (ownLabel != -1) && (neiLabel != -1) )
        {
            ++nFaces;
            ++nFacesInCell[ownLabel];
            ++nFacesInCell[neiLabel];
        }
        else if( ownLabel != -1 )
        {
            ++nFaces;
            ++nFacesInCell[ownLabel];

            if( (nei != -1) && (cType[nei] & meshOctreeAddressing::MESHCELL) )
            {
                const label procNo = octree.returnLeaf(nei).procNo();

                procFaces[procNo].append(faceI);
            }
        }
        else if( neiLabel != -1 )
        {
            ++nFaces;
            ++nFacesInCell[neiLabel];

            if( (own != -1) && (cType[own] & meshOctreeAddressing::MESHCELL) )
            {
                const label procNo = octree.returnLeaf(own).procNo();

                procFaces[procNo].append(faceI);
            }
        }
    }

    //- case is a serial run
    faces.setSize(nFaces);
    forAll(cells, cI)
        cells[cI].setSize(nFacesInCell[cI]);
    nFacesInCell = 0;

    //- calculate faces in processor patches
    if( Pstream::parRun() )
    {
        PtrList<processorBoundaryPatch>& procBoundaries =
            meshModifier.procBoundariesAccess();

        //- set the number of procBoundaries
        procBoundaries.setSize(procFaces.size());
        std::ostringstream ss;
        ss << Pstream::myProcNo();
        const word name("processor"+ss.str()+"to");
        label nProcBoundaries(nFaces), patchI(0);

        //- allocate memory for processor patches
        std::map<label, labelLongList>::const_iterator iter;
        for(iter=procFaces.begin();iter!=procFaces.end();++iter)
        {
            const label procI = iter->first;

            std::ostringstream ssNei;
            ssNei << procI;
            procBoundaries.set
            (
                patchI,
                new processorBoundaryPatch
                (
                    name+ssNei.str(),
                    "processor",
                    iter->second.size(),
                    0,
                    Pstream::myProcNo(),
                    procI
                )
            );

            nProcBoundaries -= iter->second.size();
            ++patchI;
        }

        //- create processor faces
        //- they need to be created here because of the correct ordering
        patchI = 0;
        for(iter=procFaces.begin();iter!=procFaces.end();++iter)
        {
            procBoundaries[patchI].patchStart() = nProcBoundaries;

            const labelLongList& patchFaces = iter->second;

            forAll(patchFaces, pfI)
            {
                const label fLabel = patchFaces[pfI];
                const label own = owner[fLabel];
                const label nei = neighbour[fLabel];

                const label curCell = leafCellLabel[own];
                label neiCell(-1);
                if( nei != -1 )
                    neiCell = leafCellLabel[nei];

                //- create a processor face
                if( neiCell == -1 )
                {
                    //- add a face
                    faces[nProcBoundaries].setSize(octreeFaces[fLabel].size());
                    forAllRow(octreeFaces, fLabel, pI)
                        faces[nProcBoundaries][pI] = octreeFaces(fLabel, pI);
                    cells[curCell][nFacesInCell[curCell]++] = nProcBoundaries++;
                }
                else if( curCell == -1 )
                {
                    //- add a reversed face
                    faces[nProcBoundaries].setSize(octreeFaces[fLabel].size());
                    label i(0);
                    faces[nProcBoundaries][i++] = octreeFaces(fLabel, 0);
                    for(label pI=octreeFaces.sizeOfRow(fLabel)-1;pI>0;--pI)
                        faces[nProcBoundaries][i++] = octreeFaces(fLabel, pI);
                    cells[neiCell][nFacesInCell[neiCell]++] = nProcBoundaries++;
                }
                else
                {
                    FatalErrorIn
                    (
                        "void cartesianMeshExtractor::createPolyMesh()"
                    ) << "Face " << octreeFaces[fLabel] << " causes problems!"
                        << abort(FatalError);
                }
            }

            if( procBoundaries[patchI].patchSize() !=
                (nProcBoundaries - procBoundaries[patchI].patchStart())
            )
                FatalErrorIn
                (
                    "cartesianMeshExtractor::createPolyMesh()"
                ) << "Invalid patch size!" << Pstream::myProcNo()
                    << abort(FatalError);

            ++patchI;
        }
    }

    nFaces = 0;

    forAll(octreeFaces, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        const label ownLabel = leafCellLabel[own];
        label neiLabel(-1);
        if( nei != -1 )
            neiLabel = leafCellLabel[nei];

        if( (ownLabel != -1) && (neiLabel != -1) )
        {
            //- internal face
            faces[nFaces].setSize(octreeFaces.sizeOfRow(faceI));
            forAllRow(octreeFaces, faceI, pI)
                faces[nFaces][pI] = octreeFaces(faceI, pI);

            cells[ownLabel][nFacesInCell[ownLabel]++] = nFaces;
            cells[neiLabel][nFacesInCell[neiLabel]++] = nFaces;
            ++nFaces;
        }
        else if( ownLabel != -1 )
        {
            if( (nei != -1) && (cType[nei] & meshOctreeAddressing::MESHCELL) )
            {
                //- face at a parallel boundary
                continue;
            }

            //- boundary face
            faces[nFaces].setSize(octreeFaces.sizeOfRow(faceI));
            forAllRow(octreeFaces, faceI, pI)
                faces[nFaces][pI] = octreeFaces(faceI, pI);

            cells[ownLabel][nFacesInCell[ownLabel]++] = nFaces;
            ++nFaces;
        }
        else if( neiLabel != -1 )
        {
            if( (own != -1) && (cType[own] & meshOctreeAddressing::MESHCELL) )
            {
                //- face at a parallel boundary
                continue;
            }

            //- boundary face
            faces[nFaces].setSize(octreeFaces.sizeOfRow(faceI));
            faces[nFaces][0] = octreeFaces(faceI, 0);
            for(label pI=octreeFaces.sizeOfRow(faceI)-1;pI>0;--pI)
                faces[nFaces][octreeFaces.sizeOfRow(faceI)-pI] =
                    octreeFaces(faceI, pI);

            cells[neiLabel][nFacesInCell[neiLabel]++] = nFaces;
            ++nFaces;
        }
    }

    # ifdef DEBUGMesh
    label nProcBoundaries(0);
    forAll(procBoundaries, patchI)
        nProcBoundaries += procBoundaries[patchI].patchSize();

    if( faces.size() != (nProcBoundaries + nFaces) )
    {
        Serr << "Number of faces " << faces.size() << endl;
        Serr << "Number of processor boundaries " << nProcBoundaries << endl;
        Serr << "Number of domain faces " << nFaces << endl;
        FatalErrorIn
        (
            "void cartesianMeshExtractor::createPolyMesh()"
        ) << Pstream::myProcNo() << "This mesh is invalid!"
            << abort(FatalError);
    }

    vectorField closedness(cells.size(), vector::zero);
    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();
    forAll(owner, faceI)
        if( owner[faceI] == -1 )
        {
            Info << "faces " << faces << endl;
            FatalErrorIn
            (
                "void cartesianMeshExtractor::createPolyMesh"
                "("
                "pointFieldPMG& points,"
                "faceListPMG& faces,"
                "cellListPMG& cells"
                ")"
            ) << "Face " << faceI
                << " has no owner and neighbour!!" << abort(FatalError);
        }

    forAll(faces, faceI)
    {
        const vector area = faces[faceI].normal(mesh_.points());
        closedness[owner[faceI]] += area;
        if( neighbour[faceI] != -1 )
            closedness[neighbour[faceI]] -= area;
    }

    forAll(closedness, cellI)
        if( mag(closedness[cellI]) > 1e-10 )
            Info << "Cell " << cellI << " is not closed by "
                << closedness[cellI] << endl;

    # endif

    meshModifier.reorderBoundaryFaces();

    if( octree.isQuadtree() )
    {
        //- generate empty patches
        //- search for faces with a dominant z coordinate and store them
        //- into an empty patch
        meshSurfaceEngine mse(mesh_);
        const vectorField& fNormals = mse.faceNormals();
        const faceList::subList& bFaces = mse.boundaryFaces();
        const labelList& fOwner = mse.faceOwners();
        const vectorField& fCentres = mse.faceCentres();

        const boundBox& bb = octree.rootBox();
        const scalar tZ = 0.05 * (bb.max().z() - bb.min().z());

        wordList patchNames(3);
        patchNames[0] = "defaultFaces";
        patchNames[1] = "unusedFacesBottom";
        patchNames[2] = "unusedFacesTop";

        VRWGraph boundaryFaces;
        labelLongList newFaceOwner;
        labelLongList newFacePatch;

        forAll(fNormals, bfI)
        {
            //- store the face and its owner
            boundaryFaces.appendList(bFaces[bfI]);
            newFaceOwner.append(fOwner[bfI]);

            const vector& fNormal = fNormals[bfI];

            if( Foam::mag(fNormal.z()) > Foam::mag(fNormal.x() + fNormal.y()) )
            {
                if( Foam::mag(fCentres[bfI].z() - bb.min().z()) < tZ )
                {
                    newFacePatch.append(1);
                }
                else if( Foam::mag(fCentres[bfI].z() - bb.max().z()) < tZ )
                {
                    newFacePatch.append(2);
                }
                else
                {
                    FatalErrorIn
                    (
                        "void cartesianMeshExtractor::createPolyMesh()"
                    ) << "Cannot distribute the face!!" << exit(FatalError);
                }
            }
            else
            {
                newFacePatch.append(0);
            }
        }

        //- replace the boundary with faces in correct patches
        meshModifier.replaceBoundary
        (
            patchNames,
            boundaryFaces,
            newFaceOwner,
            newFacePatch
        );

        meshModifier.boundariesAccess()[1].patchType() = "empty";
        meshModifier.boundariesAccess()[2].patchType() = "empty";
    }

    Info << "Finished creating polyMesh" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
