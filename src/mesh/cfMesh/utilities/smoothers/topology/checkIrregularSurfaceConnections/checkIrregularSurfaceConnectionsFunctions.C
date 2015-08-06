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

#include "checkIrregularSurfaceConnections.H"
#include "polyMeshGenAddressing.H"
#include "helperFunctionsPar.H"
#include "sortEdgesIntoChains.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGCheck

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool checkIrregularSurfaceConnections::checkAndFixCellGroupsAtBndVertices
(
    labelHashSet& badVertices,
    const bool removeConnections
)
{
    Info << "Checking cells connected to surface vertices" << endl;

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();

    polyMeshGenModifier meshModifier(mesh_);
    pointFieldPMG& points = meshModifier.pointsAccess();
    faceListPMG& faces = meshModifier.facesAccess();
    cellListPMG& cells = meshModifier.cellsAccess();

    const VRWGraph& cellCells = mesh_.addressingData().cellCells();
    const VRWGraph& pointCells = mesh_.addressingData().pointCells();

    boolList parallelBndNode(bPoints.size(), false);
    if( Pstream::parRun() )
    {
        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();

        forAllConstIter(Map<label>, globalToLocal, it)
            parallelBndNode[it()] = true;
    }

    //- check cells connected to surface nodes
    //- cells connected to a vertex must create a single loop when cells
    //- are visited over faces from to the other
    label nBadVertices(0);
    DynList<label> frontCells;

    const label size = bPoints.size();
    # ifdef USE_OMP
    # pragma omp parallel for private(frontCells) schedule(dynamic, 1000)
    # endif
    for(label bpI=0;bpI<size;++bpI)
    {
        if( parallelBndNode[bpI] )
            continue;

        const label pointI = bPoints[bpI];

        Map<label> cellGroup(pointCells.sizeOfRow(pointI));

        label nGroup(0);

        forAllRow(pointCells, pointI, cI)
        {
            const label cellI = pointCells(pointI, cI);

            if( cellGroup.found(cellI) )
                continue;

            cellGroup.insert(cellI, nGroup);
            frontCells.clear();
            frontCells.append(cellI);

            while( frontCells.size() != 0 )
            {
                const label cLabel = frontCells.removeLastElement();

                forAllRow(cellCells, cLabel, nI)
                {
                    const label neiCell = cellCells(cLabel, nI);

                    if( cellGroup.found(neiCell) )
                        continue;

                    if( pointCells.contains(pointI, neiCell) )
                    {
                        cellGroup.insert(neiCell, nGroup);
                        frontCells.append(neiCell);
                    }
                }

            }

            ++nGroup;
        }

        if( nGroup != 1 )
        {
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            {
                ++nBadVertices;
                badVertices.insert(pointI);

                if( removeConnections )
                {
                    const label nPoints = points.size();
                    forAllRow(pointCells, pointI, pcI)
                    {
                        const label cellI = pointCells(pointI, pcI);

                        if( cellGroup[cellI] != 0 )
                        {
                            const cell& c = cells[cellI];

                            forAll(c, fI)
                            {
                                face& f = faces[c[fI]];

                                const label pos = f.which(pointI);

                                if( pos > -1 )
                                    f[pos] = nPoints + cellGroup[cellI] - 1;
                            }
                        }
                    }

                    for(label i=1;i<nGroup;++i)
                    {
                        const point p = points[pointI];
                        points.append(p);
                    }
                }
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- check if the vertices at processor boundaries
        //- are connected correctly
        const labelList& bp = mse.bp();
        const label origNumVertices = bp.size();

        const labelList& owner = mesh_.owner();
        const labelList& neighbour = mesh_.neighbour();
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh_.procBoundaries();

        const labelLongList& globalCellLabel =
            mesh_.addressingData().globalCellLabel();

        std::map<label, DynList<edge> > dualEdgesForPoint;
        std::map<label, DynList<edge> >::iterator bpIter;
        forAll(parallelBndNode, bpI)
        {
            if( !parallelBndNode[bpI] )
                continue;

            dualEdgesForPoint.insert
            (
                std::make_pair(bpI, DynList<edge>())
            );
        }

        //- fill-in dualEdgesForPoint with local data
        for(label faceI=0;faceI<mesh_.nInternalFaces();++faceI)
        {
            const face& f = faces[faceI];

            forAll(f, pI)
            {
                if( f[pI] >= origNumVertices )
                    continue;

                bpIter = dualEdgesForPoint.find(bp[f[pI]]);
                if( bpIter != dualEdgesForPoint.end() )
                {
                    const label cOwn = globalCellLabel[owner[faceI]];
                    const label cNei = globalCellLabel[neighbour[faceI]];

                    //- store the edge
                    bpIter->second.append(edge(cOwn, cNei));
                }
            }
        }

        //- fill-in with data at processor boundaries. Store edges
        //- on the processor with the lower label not to duplicate the data
        forAll(procBoundaries, patchI)
        {
            if( procBoundaries[patchI].owner() )
                continue;

            const label start = procBoundaries[patchI].patchStart();
            labelList globalLabels(procBoundaries[patchI].patchSize());
            forAll(globalLabels, fI)
                globalLabels[fI] = globalCellLabel[owner[start+fI]];

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                globalLabels.byteSize()
            );

            toOtherProc << globalLabels;
        }

        forAll(procBoundaries, patchI)
        {
            if( !procBoundaries[patchI].owner() )
                continue;

            labelList receivedData;
            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            fromOtherProc >> receivedData;

            const label start = procBoundaries[patchI].patchStart();

            forAll(receivedData, i)
            {
                const face& f = faces[start+i];

                forAll(f, pI)
                {
                    bpIter = dualEdgesForPoint.find(bp[f[pI]]);
                    if( bpIter != dualEdgesForPoint.end() )
                    {
                        const label cOwn = globalCellLabel[owner[start+i]];
                        const label cNei = receivedData[i];
                        bpIter->second.append(edge(cOwn, cNei));
                    }
                }
            }
        }

        //- exchange data with other processors
        //- this step supplies all processors with all necessary data
        const VRWGraph& bpAtProcs = mse.bpAtProcs();
        const labelList& globalPointLabel = mse.globalBoundaryPointLabel();
        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();
        std::map<label, labelLongList> exchangeData;
        for
        (
            bpIter=dualEdgesForPoint.begin();
            bpIter!=dualEdgesForPoint.end();
            ++bpIter
        )
        {
            const label bpI = bpIter->first;

            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);
                if( neiProc == Pstream::myProcNo() )
                    continue;

                //- edges are sent to other processors as follows
                //- 1. global point label
                //- 2. number of edges at node
                //- 3. labels of edges
                labelLongList& dts = exchangeData[neiProc];
                const DynList<edge>& edges = bpIter->second;
                dts.append(globalPointLabel[bpI]);
                dts.append(edges.size());
                forAll(edges, eI)
                {
                    dts.append(edges[eI].start());
                    dts.append(edges[eI].end());
                }
            }
        }

        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        label counter = 0;
        while( counter < receivedData.size() )
        {
            const label bpI = globalToLocal[receivedData[counter++]];
            const label nEdges = receivedData[counter++];

            for(label i=0;i<nEdges;++i)
            {
                const label s = receivedData[counter++];
                const label e = receivedData[counter++];
                dualEdgesForPoint[bpI].append(edge(s, e));
            }
        }

        # ifdef DEBUGCheck
        for(label i=0;i<Pstream::nProcs();++i)
        {
            if( i == Pstream::myProcNo() )
            {
                for
                (
                    bpIter=dualEdgesForPoint.begin();
                    bpIter!=dualEdgesForPoint.end();
                    ++bpIter
                )
                {
                    const DynList<edge>& pEdges = bpIter->second;
                    Pout << "Global point label " << globalPointLabel[bpIter->first]
                        << " point cells " << pointCells[bPoints[bpIter->first]]
                        << " dual edges " << pEdges << endl;
                }
            }

            returnReduce(1, sumOp<label>());
        }
        # endif

        //- Finally, check the number of dual loops processor vertices
        for
        (
            bpIter=dualEdgesForPoint.begin();
            bpIter!=dualEdgesForPoint.end();
            ++bpIter
        )
        {
            const DynList<edge>& pEdges = bpIter->second;
            std::map<label, DynList<label> > bpEdges;
            forAll(pEdges, eI)
            {
                forAll(pEdges[eI], i)
                    bpEdges[pEdges[eI][i]].append(eI);
            }

            //- check if all points can be visited via edges
            counter = 0;
            Map<label> cellGroup(pEdges.size());
            for
            (
                std::map<label, DynList<label> >::iterator it=bpEdges.begin();
                it!=bpEdges.end();
                ++it
            )
            {
                if( cellGroup.found(it->first) )
                    continue;

                cellGroup.insert(it->first, counter);
                frontCells.clear();
                frontCells.append(it->first);

                while( frontCells.size() != 0 )
                {
                    const label cLabel = frontCells.removeLastElement();

                    const DynList<label>& attachedEdges = bpEdges[cLabel];
                    forAll(attachedEdges, aeI)
                    {
                        const label oCell =
                            pEdges[attachedEdges[aeI]].otherVertex(cLabel);

                        if( cellGroup.found(oCell) )
                            continue;

                        frontCells.append(oCell);
                        cellGroup.insert(oCell, counter);
                    }
                }

                ++counter;
            }

            if( counter != 1 )
            {
                ++nBadVertices;
                badVertices.insert(bPoints[bpIter->first]);

                if( !removeConnections )
                    continue;

                const label pointI = bPoints[bpIter->first];
                const label nPoints = points.size();
                forAllRow(pointCells, pointI, pcI)
                {
                    const label cellI = pointCells(pointI, pcI);

                    if( cellGroup[globalCellLabel[cellI]] == 0 )
                        continue;

                    const cell& c = cells[cellI];

                    forAll(c, fI)
                    {
                        face& f = faces[c[fI]];

                        const label pos = f.which(pointI);

                        if( pos > -1 )
                            f[pos] =
                                nPoints + cellGroup[globalCellLabel[cellI]] - 1;
                    }
                }

                for(label i=1;i<counter;++i)
                {
                    const point p = points[pointI];
                    points.append(p);
                }
            }
        }
    }

    reduce(nBadVertices, sumOp<label>());

    Info << "Found " << nBadVertices << " problematic vertices" << endl;
    Info << "Finished checking cells connected to surface vertices" << endl;

    if( nBadVertices != 0 )
    {
        clearMeshEngine();
        mesh_.clearAddressingData();

        return true;
    }

    return false;
}

bool checkIrregularSurfaceConnections::checkEdgeFaceConnections
(
    labelHashSet& badVertices,
    const bool removeCells
)
{
    Info << "Checking for non-manifold surface edges" << endl;

    const meshSurfaceEngine& mse = surfaceEngine();
    const labelList& bp = mse.bp();
    const edgeList& edges = mse.edges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const VRWGraph& pointEdges = mse.boundaryPointEdges();

    labelHashSet badEdges;

    forAll(edgeFaces, edgeI)
        if( edgeFaces.sizeOfRow(edgeI) > 2 )
        {
            badVertices.insert(edges[edgeI].start());
            badVertices.insert(edges[edgeI].end());

            badEdges.insert(edgeI);
        }

    if( Pstream::parRun() )
    {
        //- boundary edges at processor boundaries
        Map<label> numFacesAtEdge;
        const labelList& globalEdgeLabel = mse.globalBoundaryEdgeLabel();
        const Map<label>& globalToLocalEdgeLabel =
            mse.globalToLocalBndEdgeAddressing();
        const VRWGraph& edgesAtProcs = mse.beAtProcs();

        const DynList<label>& neiProcs = mse.beNeiProcs();
        std::map<label, labelLongList> exchangeData;
        forAll(neiProcs, procI)
            exchangeData.insert
            (
                std::make_pair(neiProcs[procI], labelLongList())
            );
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
                    eIter->second.append(globalEdgeLabel[eI]);
                    eIter->second.append(edgeFaces.sizeOfRow(eI));
                }
            }
        }

        //- send data to other processors
        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        label counter(0);
        while( counter < receivedData.size() )
        {
            const label geI = receivedData[counter++];
            const label nFaces = receivedData[counter++];

            numFacesAtEdge[geI] += nFaces;

            if( numFacesAtEdge[geI] > 2 )
            {
                const label edgeI = globalToLocalEdgeLabel[geI];
                badVertices.insert(edges[edgeI].start());
                badVertices.insert(edges[edgeI].end());

                badEdges.insert(edgeI);
            }
        }
    }

    const label nBadEdges = returnReduce(badEdges.size(), sumOp<label>());
    Info << "Found " << nBadEdges << " non-manifold edges" << endl;
    Info << "Finished checking for non-manifold surface edges" << endl;

    if( nBadEdges != 0 && removeCells )
    {
        //- remove all cells connected to the selected edge
        boolList removeCell(mesh_.cells().size(), false);

        const faceListPMG& faces = mesh_.faces();
        const labelList& owner = mesh_.owner();
        const labelList& neighbour = mesh_.neighbour();

        forAll(faces, faceI)
        {
            const face& f = faces[faceI];

            forAll(f, pI)
            {
                const edge e = f.faceEdge(pI);

                if( (bp[e[0]] != -1) && (bp[e[1]] != -1) )
                {
                    const label bpI = bp[e[0]];

                    forAllRow(pointEdges, bpI, peI)
                    {
                        const label edgeI = pointEdges(bpI, peI);
                        if( (edges[edgeI] == e) && badEdges.found(edgeI) )
                        {
                            removeCell[owner[faceI]] = true;
                            if( neighbour[faceI] < 0 )
                                continue;
                            removeCell[neighbour[faceI]] = true;
                        }
                    }
                }
            }
        }

        polyMeshGenModifier(mesh_).removeCells(removeCell);
        clearMeshEngine();

        return true;
    }

    return false;
}

bool checkIrregularSurfaceConnections::checkFaceGroupsAtBndVertices
(
    labelHashSet& badVertices,
    const bool removeCells
)
{
    Info << "Checking faces connections to surface vertices" << endl;

    labelHashSet invalidVertices(100);

    const meshSurfaceEngine& mse = surfaceEngine();

    const labelList& bPoints = mse.boundaryPoints();
    const VRWGraph& pointFaces = mse.pointFaces();
    const VRWGraph& faceFaces = mse.faceFaces();

    boolList parallelBndPoint(bPoints.size(), false);
    if( Pstream::parRun() )
    {
        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();

        forAllConstIter(Map<label>, globalToLocal, it)
            parallelBndPoint[it()] = true;
    }

    //- check number of face groups
    DynList<label> front;
    const label size = pointFaces.size();
    # ifdef USE_OMP
    # pragma omp parallel for private(front) schedule(dynamic)
    # endif
    for(label bpI=0;bpI<size;++bpI)
    {
        if( parallelBndPoint[bpI] )
            continue;

        Map<label> faceGroup(pointFaces.sizeOfRow(bpI));
        label nGroup(0);

        forAllRow(pointFaces, bpI, pfI)
        {
            const label fI = pointFaces(bpI, pfI);

            if( faceGroup.found(fI) )
                continue;

            front.clear();
            front.append(fI);
            faceGroup.insert(fI, nGroup);

            while( front.size() != 0 )
            {
                const label fLabel = front.removeLastElement();

                forAllRow(faceFaces, fLabel, ffI)
                {
                    const label neiFace = faceFaces(fLabel, ffI);

                    if( faceGroup.found(neiFace) )
                        continue;

                    if( pointFaces.contains(bpI, neiFace) )
                    {
                        front.append(neiFace);
                        faceGroup.insert(neiFace, nGroup);
                    }
                }
            }

            ++nGroup;
        }

        if( nGroup != 1 )
        {
            const label pointI = bPoints[bpI];
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            {
                badVertices.insert(pointI);
                invalidVertices.insert(pointI);
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- check connections at parallel vertices
        //- a connection of two faces over an edge can be represented
        //- as an edge. A list of edges at a bnd vertex must be connected
        //- into a single loop, otherwise the surface is ill-connected.

        const labelList& globalPointLabel = mse.globalBoundaryPointLabel();
        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();
        const VRWGraph& bpAtProcs = mse.bpAtProcs();
        const labelList& bp = mse.bp();

        const edgeList& edges = mse.edges();
        const VRWGraph& pointEdges = mse.boundaryPointEdges();
        const VRWGraph& edgeFaces = mse.edgeFaces();

        const labelList& globalFaceLabel = mse.globalBoundaryFaceLabel();

        const labelList& globalEdgeLabel = mse.globalBoundaryEdgeLabel();
        const Map<label>& globalToLocalEdge =
            mse.globalToLocalBndEdgeAddressing();
        const Map<label>& otherFaceAtProc = mse.otherEdgeFaceAtProc();
        const DynList<label>& beNeiProcs = mse.beNeiProcs();

        //- create map of dual edges for boundary points at processor boundaries
        std::map<label, DynList<edge> > dualEdgesForPoint;
        std::map<label, DynList<edge> >::iterator bpIter;
        forAll(parallelBndPoint, bpI)
        {
            if( !parallelBndPoint[bpI] )
                continue;

            dualEdgesForPoint.insert(std::make_pair(bpI, DynList<edge>()));

            forAllRow(pointEdges, bpI, peI)
            {
                const label edgeI = pointEdges(bpI, peI);

                if( edgeFaces.sizeOfRow(edgeI) == 2 )
                {
                    //- add local edge
                    edge e
                    (
                        globalFaceLabel[edgeFaces(edgeI, 0)],
                        globalFaceLabel[edgeFaces(edgeI, 1)]
                    );

                    dualEdgesForPoint[bpI].append(e);
                }
            }
        }

        std::map<label, labelLongList> exchangeData;
        //- collect connections over processor edges on the processor with
        //- the lowest label to avoid duplication of data
        forAll(beNeiProcs, i)
            exchangeData.insert(std::make_pair(beNeiProcs[i], labelLongList()));
        forAllConstIter(Map<label>, otherFaceAtProc, it)
        {
            const label beI = it.key();

            if( it() < Pstream::myProcNo() )
            {
                //- the data is sent as follows:
                //- 1. global edge label
                //- 2. global face label
                exchangeData[it()].append(globalEdgeLabel[beI]);
                exchangeData[it()].append(globalFaceLabel[edgeFaces(beI, 0)]);
            }
        }

        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        label counter(0);
        while( counter < receivedData.size() )
        {
            const label beI = globalToLocalEdge[receivedData[counter++]];
            const label fLabel = receivedData[counter++];

            const edge& e = edges[beI];
            forAll(e, i)
            {
                const label bpI = bp[e[i]];

                const label fLocal = globalFaceLabel[edgeFaces(beI, 0)];
                dualEdgesForPoint[bpI].append(edge(fLabel, fLocal));
            }
        }

        //- exchange data with other processors
        //- this step supplies all processors with all necessary data
        exchangeData.clear();
        for
        (
            bpIter=dualEdgesForPoint.begin();
            bpIter!=dualEdgesForPoint.end();
            ++bpIter
        )
        {
            const label bpI = bpIter->first;

            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);
                if( neiProc == Pstream::myProcNo() )
                    continue;

                //- edges are sent to other processors as follows
                //- 1. global point label
                //- 2. number of edges at node
                //- 3. labels of edges
                labelLongList& dts = exchangeData[neiProc];
                const DynList<edge>& edges = bpIter->second;
                dts.append(globalPointLabel[bpI]);
                dts.append(edges.size());
                forAll(edges, eI)
                {
                    dts.append(edges[eI].start());
                    dts.append(edges[eI].end());
                }
            }
        }

        receivedData.clear();
        help::exchangeMap(exchangeData, receivedData);

        counter = 0;
        while( counter < receivedData.size() )
        {
            const label bpI = globalToLocal[receivedData[counter++]];
            const label nEdges = receivedData[counter++];
            for(label i=0;i<nEdges;++i)
            {
                const label s = receivedData[counter++];
                const label e = receivedData[counter++];
                dualEdgesForPoint[bpI].append(edge(s, e));
            }
        }

        # ifdef DEBUGCheck
        for(label i=0;i<Pstream::nProcs();++i)
        {
            if( i == Pstream::myProcNo() )
            {
                for
                (
                    bpIter=dualEdgesForPoint.begin();
                    bpIter!=dualEdgesForPoint.end();
                    ++bpIter
                )
                {
                    const DynList<edge>& pEdges = bpIter->second;
                    Pout << "Global point label "
                        << globalPointLabel[bpIter->first]
                        << " dual edges " << pEdges << endl;
                }
            }

            returnReduce(1, sumOp<label>());
        }
        # endif

        //- Finally, check the number of dual loops processor vertices
        for
        (
            bpIter=dualEdgesForPoint.begin();
            bpIter!=dualEdgesForPoint.end();
            ++bpIter
        )
        {
            const DynList<edge>& pEdges = bpIter->second;

            std::map<label, DynList<label> > bpEdges;
            forAll(pEdges, eI)
            {
                forAll(pEdges[eI], i)
                    bpEdges[pEdges[eI][i]].append(eI);
            }

            //- check if all points can be visited via edges
            counter = 0;
            Map<label> cellGroup(pEdges.size());
            for
            (
                std::map<label, DynList<label> >::iterator it=bpEdges.begin();
                it!=bpEdges.end();
                ++it
            )
            {
                if( cellGroup.found(it->first) )
                    continue;

                cellGroup.insert(it->first, counter);
                front.clear();
                front.append(it->first);

                while( front.size() != 0 )
                {
                    const label cLabel = front.removeLastElement();

                    const DynList<label>& attachedEdges = bpEdges[cLabel];
                    forAll(attachedEdges, aeI)
                    {
                        const label oCell =
                            pEdges[attachedEdges[aeI]].otherVertex(cLabel);

                        if( cellGroup.found(oCell) )
                            continue;

                        front.append(oCell);
                        cellGroup.insert(oCell, counter);
                    }
                }

                ++counter;
            }

            if( counter != 1 )
            {
                invalidVertices.insert(bPoints[bpIter->first]);
                badVertices.insert(bPoints[bpIter->first]);
            }
        }
    }

    const label nBadVertices =
        returnReduce(invalidVertices.size(), sumOp<label>());
    Info << "Found " << nBadVertices << " invalid connections" << endl;
    Info << "Finished checking faces connections to surface vertices" << endl;

    if( nBadVertices != 0 && removeCells )
    {
        const VRWGraph& pointCells = mesh_.addressingData().pointCells();

        boolList removeCell(mesh_.cells().size(), false);

        forAllConstIter(labelHashSet, invalidVertices, it)
        {
            forAllRow(pointCells, it.key(), pcI)
                removeCell[pointCells(it.key(), pcI)] = true;
        }

        polyMeshGenModifier(mesh_).removeCells(removeCell);
        clearMeshEngine();

        return true;
    }

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
