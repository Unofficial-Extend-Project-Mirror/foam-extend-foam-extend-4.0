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

#include "topologicalCleaner.H"
#include "polyMeshGenAddressing.H"
#include "DynList.H"
#include "meshSurfaceEngine.H"

#include <map>

//#define DEBUGCheck

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void topologicalCleaner::checkInvalidConnectionsForVerticesCells
(
    labelHashSet* irregularNodesPtr
)
{
    if( Pstream::parRun() )
    {
        return;

        FatalErrorIn
        (
            "void topologicalCleaner::checkInvalidConnectionsForVerticesCells()"
        ) << "This does not run in parallel!" << exit(FatalError);
    }

    polyMeshGenModifier meshModifier(mesh_);

    label nPoints = mesh_.points().size();
    pointFieldPMG& points = meshModifier.pointsAccess();
    faceListPMG& faces = meshModifier.facesAccess();
    const cellListPMG& cells = mesh_.cells();

    const VRWGraph& pointCells = mesh_.addressingData().pointCells();
    const VRWGraph& cellCells = mesh_.addressingData().cellCells();

    meshSurfaceEngine mse(mesh_);
    const labelList& bPoints = mse.boundaryPoints();

    label nInvalidConnections(0);

    forAll(bPoints, bpI)
    {
        const label pointI = bPoints[bpI];

        # ifdef DEBUGCheck
        Info << "Checking point " << pointI << endl;
        # endif

        label material(1);

        labelList materialForCell(pointCells.sizeOfRow(pointI), 0);

        forAllRow(pointCells, pointI, cI)
        if( !materialForCell[cI] )
        {
            materialForCell[cI] = material;

            DynList<label> frontCells;
            frontCells.append(cI);

            do
            {
                DynList<label> newFrontCells;

                forAll(frontCells, fcI)
                {
                    const label pointCellI =
                    pointCells(pointI, frontCells[fcI]);

                    forAllRow(cellCells, pointCellI, nI)
                    {
                        forAllRow(pointCells, pointI, pcI)
                        {
                            if( materialForCell[pcI] )
                                continue;

                            if(
                                cellCells(pointCellI, nI) ==
                                pointCells(pointI, pcI)
                            )
                            {
                                newFrontCells.append(pcI);
                                materialForCell[pcI] = material;
                                break;
                            }
                        }
                    }
                }

                frontCells = newFrontCells;

            } while( frontCells.size() != 0 );

            ++material;
        }

        # ifdef DEBUGCheck
        Info << "Number of materials for vertex is " << material << endl;
        # endif

        if( material > 2 )
        {
            ++nInvalidConnections;

            if( irregularNodesPtr )
            {
                irregularNodesPtr->insert(pointI);
                continue;
            }

            forAllRow(pointCells, pointI, pcI)
            {
                if( materialForCell[pcI] == 1 )
                    continue;

                const cell& c = cells[pointCells(pointI, pcI)];

                forAll(c, fI)
                {
                    face& f = faces[c[fI]];

                    forAll(f, pI)
                        if( f[pI] == pointI )
                            f[pI] = nPoints + materialForCell[pcI] - 1;
                }
            }

            for(label i=1;i<material;++i)
            {
                const point p = points[pointI];
                points.append(p);
                ++nPoints;
            }
        }
    }

    Info << "Found " << nInvalidConnections
        << " invalid cell connections" << endl;

    mesh_.clearAddressingData();

    if( nInvalidConnections != 0 )
        meshModifier.removeUnusedVertices();
}

void topologicalCleaner::checkInvalidConnectionsForVerticesFaces
(
    labelHashSet* irregularNodesPtr
)
{
    const meshSurfaceEngine mse(mesh_);

    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelList& faceOwner = mse.faceOwners();

    boolList removeCell(mesh_.cells().size(), false);
    bool changed(false);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(static, 1)
    # endif
    forAll(edgeFaces, edgeI)
        if( edgeFaces.sizeOfRow(edgeI) > 2 )
        {
            forAllRow(edgeFaces, edgeI, fI)
                removeCell[faceOwner[edgeFaces(edgeI, fI)]] = true;

            changed = true;
        }

    if( Pstream::parRun() )
    {
        //- boundary edges at processor boundaries
        Map<label> numFacesAtEdge;
        const labelList& globalEdgeLabel = mse.globalBoundaryEdgeLabel();
        const Map<label>& globalToLocal = mse.globalToLocalBndEdgeAddressing();
        const VRWGraph& edgesAtProcs = mse.beAtProcs();

        DynList<label> neiProcs;
        std::map<label, labelLongList> exchangeData;
        std::map<label, labelLongList>::iterator eIter;

        forAll(edgeFaces, eI)
        {
            if( edgesAtProcs.sizeOfRow(eI) > 0 )
            {
                numFacesAtEdge.insert
                (
                    globalEdgeLabel[eI],
                    edgeFaces.sizeOfRow(eI)
                );

                forAllRow(edgesAtProcs, eI, procI)
                {
                    const label neiProc = edgesAtProcs(eI, procI);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    eIter = exchangeData.find(neiProc);
                    if( eIter == exchangeData.end() )
                    {
                        neiProcs.append(neiProc);
                        exchangeData.insert
                        (
                            std::pair<label, labelLongList>
                            (
                                neiProc,
                                labelLongList()
                            )
                        );

                        eIter = exchangeData.find(neiProc);
                    }

                    eIter->second.append(globalEdgeLabel[eI]);
                    eIter->second.append(edgeFaces.sizeOfRow(eI));
                }
            }
        }

        //- send data to other processors
        forAll(neiProcs, procI)
        {
            eIter = exchangeData.find(neiProcs[procI]);
            const labelLongList& dataToSend = eIter->second;

            OPstream toOtherProc
            (
                Pstream::blocking,
                neiProcs[procI],
                dataToSend.byteSize()
            );
            toOtherProc << dataToSend;
        }

        forAll(neiProcs, procI)
        {
            labelList receivedData;
            IPstream fromOtherProc(Pstream::blocking, neiProcs[procI]);
                fromOtherProc >> receivedData;

            label counter(0);
            while( counter < receivedData.size() )
            {
                const label geI = receivedData[counter++];
                const label nFaces = receivedData[counter++];

                numFacesAtEdge[geI] += nFaces;

                if( numFacesAtEdge[geI] > 2 )
                {
                    const label edgeI = globalToLocal[geI];
                    forAllRow(edgeFaces, edgeI, fI)
                        removeCell[faceOwner[edgeFaces(edgeI, fI)]] = true;

                    changed = true;
                }
            }
        }
    }

    reduce(changed, maxOp<bool>());

    if( changed )
    {
        polyMeshGenModifier(mesh_).removeCells(removeCell);

        decomposeCell_.setSize(mesh_.cells().size());
        decomposeCell_ = false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
