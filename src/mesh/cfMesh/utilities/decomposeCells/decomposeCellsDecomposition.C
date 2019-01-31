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

#include "decomposeCells.H"
#include "demandDrivenData.H"
#include "polyMeshGenAddressing.H"
#include "meshSurfaceEngine.H"
#include "decomposeFaces.H"
#include "labelLongList.H"

//#define DEBUGDecompose

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void decomposeCells::decomposeMesh(const boolList& decomposeCell)
{
    checkFaceConnections(decomposeCell);

    createPointsAndCellFaces(decomposeCell);

    storeBoundaryFaces(decomposeCell);

    removeDecomposedCells(decomposeCell);

    addNewCells();

    # ifdef DEBUGDecompose
    mesh_.addressingData().checkMesh();
    # endif
}

void decomposeCells::checkFaceConnections(const boolList& decomposeCell)
{
    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();

    boolList decomposeFace(faces.size(), false);
    forAll(cells, cellI)
    {
        if( decomposeCell[cellI] )
        {
            DynList<label, 32> vrt;
            DynList<edge, 64> edges;
            DynList<DynList<label, 8> > faceEdges;
            DynList<DynList<label, 2>, 64> edgeFaces;

            findAddressingForCell(cellI, vrt, edges, faceEdges, edgeFaces);

            forAll(faceEdges, fI)
            {
                const DynList<label, 8>& fEdges = faceEdges[fI];

                labelHashSet neiFaces;
                forAll(fEdges, feI)
                {
                    label neiFace = edgeFaces[fEdges[feI]][0];
                    if( neiFace == fI )
                        neiFace = edgeFaces[fEdges[feI]][1];

                    if( neiFaces.found(neiFace) )
                    {
                        decomposeFace[cells[cellI][fI]] = true;
                    }
                    else
                    {
                        neiFaces.insert(neiFace);
                    }
                }
            }
        }
    }

    if( Pstream::parRun() )
    {
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh_.procBoundaries();

        //- send information to the neighbour processor
        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            boolList decFace(procBoundaries[patchI].patchSize(), false);
            const label size = decFace.size();

            for(label i=0;i<size;++i)
            {
                if( decomposeFace[start+i] )
                    decFace[i] = true;
            }

            OPstream toOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo(),
                decFace.byteSize()
            );

            toOtherProc << decFace;
        }

        //- receive information from the neighbour processor
        forAll(procBoundaries, patchI)
        {
            boolList decFace;

            IPstream fromOtherProc
            (
                Pstream::commsTypes::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            fromOtherProc >> decFace;

            const label start = procBoundaries[patchI].patchStart();
            forAll(decFace, i)
            {
                if( decFace[i] )
                    decomposeFace[start+i] = true;
            }
        }
    }

    //- decompose faces which would cause invalid connections
    decomposeFaces(mesh_).decomposeMeshFaces(decomposeFace);
}

void decomposeCells::createPointsAndCellFaces(const boolList& decomposeCell)
{
    facesOfNewCells_.clear();

    forAll(decomposeCell, cI)
        if( decomposeCell[cI] )
        {
            decomposeCellIntoPyramids(cI);
        }
}

void decomposeCells::storeBoundaryFaces(const boolList& /*decomposeCell*/)
{
    meshSurfaceEngine mse(mesh_);
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& facePatch = mse.boundaryFacePatches();

    forAll(bFaces, bfI)
    {
        newBoundaryFaces_.appendList(bFaces[bfI]);
        newBoundaryPatches_.append(facePatch[bfI]);
    }
}

void decomposeCells::removeDecomposedCells(const boolList& decomposeCell)
{
    # ifdef DEBUGDecompose
    Info << "Number of cells before removal " << mesh_.cells().size() << endl;
    # endif

    polyMeshGenModifier meshModifier(mesh_);
    meshModifier.removeCells(decomposeCell, false);

    # ifdef DEBUGDecompose
    Info << "Number of cells after removal " << mesh_.cells().size() << endl;
    # endif
}

void decomposeCells::addNewCells()
{
    Info << "Adding new cells " << endl;
    polyMeshGenModifier(mesh_).addCells(facesOfNewCells_);
    facesOfNewCells_.clear();
    Info << "Reordering bnd faces" << endl;
    polyMeshGenModifier(mesh_).reorderBoundaryFaces();

    Info << "Finding bnd faces" << endl;
    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();
    const VRWGraph& pointFaces = mesh_.addressingData().pointFaces();

    labelLongList newBoundaryOwners;

    forAll(newBoundaryFaces_, faceI)
    {
        face bf(newBoundaryFaces_.sizeOfRow(faceI));
        forAllRow(newBoundaryFaces_, faceI, pI)
            bf[pI] = newBoundaryFaces_(faceI, pI);

        # ifdef DEBUGDecompose
        Info << "Finding cell for boundary face " << bf << endl;
        bool found(false);
        forAllRow(pointFaces, bf[0], pfI)
            if( bf == faces[pointFaces(bf[0], pfI)] )
                found = true;
        if( !found )
            FatalErrorIn
            (
                "void decomposeCells::addNewCells()"
            ) << "Face " << bf << " does not exist in the mesh"
                << abort(FatalError);
        #endif

        forAllRow(pointFaces, bf[0], pfI)
        {
            const label fLabel = pointFaces(bf[0], pfI);
            if( (mesh_.faceIsInPatch(fLabel) != -1) && (bf == faces[fLabel]) )
            {
                # ifdef DEBUGDecompose
                Info << "Boundary face " << bf << " is in cell "
                << owner[fLabel]] << endl;
                # endif

                newBoundaryOwners.append(owner[fLabel]);
            }
        }
    }

    polyMeshGenModifier(mesh_).replaceBoundary
    (
        patchNames_,
        newBoundaryFaces_,
        newBoundaryOwners,
        newBoundaryPatches_
    );

    polyMeshGenModifier(mesh_).removeUnusedVertices();
    polyMeshGenModifier(mesh_).clearAll();

    PtrList<boundaryPatch>& boundaries =
        polyMeshGenModifier(mesh_).boundariesAccess();
    forAll(boundaries, patchI)
        boundaries[patchI].patchType() = patchTypes_[patchI];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// ************************************************************************* //
