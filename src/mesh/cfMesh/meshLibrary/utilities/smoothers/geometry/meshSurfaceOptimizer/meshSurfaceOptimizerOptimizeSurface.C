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

#include "demandDrivenData.H"
#include "meshSurfaceOptimizer.H"
#include "meshSurfaceEngineModifier.H"
#include "meshSurfaceCheckInvertedVertices.H"
#include "meshOctree.H"
#include "triangle.H"
#include "helperFunctionsPar.H"
#include "meshSurfaceMapper.H"
#include "meshSurfaceMapper2D.H"
#include "polyMeshGen2DEngine.H"
#include "polyMeshGenAddressing.H"
#include "labelledPoint.H"
#include "FIFOStack.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label meshSurfaceOptimizer::findInvertedVertices
(
    boolList& smoothVertex,
    const label nAdditionalLayers
) const
{
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const VRWGraph& pPoints = surfaceEngine_.pointPoints();

    if( smoothVertex.size() != bPoints.size() )
    {
        smoothVertex.setSize(bPoints.size());
        smoothVertex = true;
    }

    label nInvertedTria(0);

    //- check the vertices at the surface
    //- mark the ones where the mesh is tangled
    meshSurfaceCheckInvertedVertices vrtCheck(surfaceEngine_, smoothVertex);
    const labelHashSet& inverted = vrtCheck.invertedVertices();

    smoothVertex = false;
    forAll(bPoints, bpI)
    {
        if( inverted.found(bPoints[bpI]) )
        {
            ++nInvertedTria;
            smoothVertex[bpI] = true;
        }
    }

    if( Pstream::parRun() )
        reduce(nInvertedTria, sumOp<label>());
    Info << "Number of inverted boundary faces is " << nInvertedTria << endl;

    if( nInvertedTria == 0 )
        return 0;

    //- add additional layers around inverted points
    for(label i=0;i<nAdditionalLayers;++i)
    {
        boolList originallySelected = smoothVertex;
        forAll(smoothVertex, bpI)
            if( originallySelected[bpI] )
                forAllRow(pPoints, bpI, ppI)
                    smoothVertex[pPoints(bpI, ppI)] = true;

        if( Pstream::parRun() )
        {
            //- exchange global labels of inverted points
            const labelList& globalPointLabel =
                surfaceEngine_.globalBoundaryPointLabel();
            const Map<label>& globalToLocal =
                surfaceEngine_.globalToLocalBndPointAddressing();
            const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
            const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();

            std::map<label, labelLongList> shareData;
            forAll(neiProcs, procI)
                shareData.insert
                (
                    std::make_pair(neiProcs[procI], labelLongList())
                );

            forAllConstIter(Map<label>, globalToLocal, iter)
            {
                const label bpI = iter();

                if( !smoothVertex[bpI] )
                    continue;

                forAllRow(bpAtProcs, bpI, procI)
                {
                    const label neiProc = bpAtProcs(bpI, procI);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    shareData[neiProc].append(globalPointLabel[bpI]);
                }
            }

            //- exchange data with other processors
            labelLongList receivedData;
            help::exchangeMap(shareData, receivedData);

            forAll(receivedData, j)
            {
                const label bpI = globalToLocal[receivedData[j]];

                smoothVertex[bpI] = true;
            }
        }
    }

    return nInvertedTria;
}

void meshSurfaceOptimizer::smoothEdgePoints
(
    const labelLongList& edgePoints,
    const labelLongList& procEdgePoints
)
{
    DynList<LongList<labelledPoint> > newPositions(1);
    # ifdef USE_OMP
    newPositions.setSize(omp_get_num_procs());
    # endif

    //- smooth edge vertices
    # ifdef USE_OMP
    # pragma omp parallel num_threads(newPositions.size())
    # endif
    {
        # ifdef USE_OMP
        LongList<labelledPoint>& newPos =
            newPositions[omp_get_thread_num()];
        # else
        LongList<labelledPoint>& newPos = newPositions[0];
        # endif

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 40)
        # endif
        forAll(edgePoints, i)
        {
            const label bpI = edgePoints[i];

            if( vertexType_[bpI] & PROCBND )
                continue;

            newPos.append(labelledPoint(bpI, newEdgePositionLaplacian(bpI)));
        }
    }

    if( Pstream::parRun() )
        edgeNodeDisplacementParallel(procEdgePoints);

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);
    forAll(newPositions, threadI)
    {
        const LongList<labelledPoint>& newPos = newPositions[threadI];

        forAll(newPos, i)
            surfaceModifier.moveBoundaryVertexNoUpdate
            (
                newPos[i].pointLabel(),
                newPos[i].coordinates()
            );
    }

    surfaceModifier.updateGeometry(edgePoints);
}

void meshSurfaceOptimizer::smoothLaplacianFC
(
    const labelLongList& selectedPoints,
    const labelLongList& selectedProcPoints,
    const bool transform
)
{
    DynList<LongList<labelledPoint> > newPositions(1);
    # ifdef USE_OMP
    newPositions.setSize(omp_get_num_procs());
    # endif

    # ifdef USE_OMP
    # pragma omp parallel num_threads(newPositions.size())
    # endif
    {
        # ifdef USE_OMP
        LongList<labelledPoint>& newPos =
            newPositions[omp_get_thread_num()];
        # else
        LongList<labelledPoint>& newPos = newPositions[0];
        # endif

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 40)
        # endif
        forAll(selectedPoints, i)
        {
            const label bpI = selectedPoints[i];

            if( vertexType_[bpI] & PROCBND )
                continue;

            newPos.append
            (
                labelledPoint(bpI, newPositionLaplacianFC(bpI, transform))
            );
        }
    }

    if( Pstream::parRun() )
        nodeDisplacementLaplacianFCParallel(selectedProcPoints, transform);

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);

    # ifdef USE_OMP
    # pragma omp parallel num_threads(newPositions.size())
    # endif
    {
        # ifdef USE_OMP
        const label threadI = omp_get_thread_num();
        # else
        const label threadI = 0;
        # endif

        const LongList<labelledPoint>& newPos = newPositions[threadI];

        forAll(newPos, i)
            surfaceModifier.moveBoundaryVertexNoUpdate
            (
                newPos[i].pointLabel(),
                newPos[i].coordinates()
            );
    }

    surfaceModifier.updateGeometry(selectedPoints);
}

void meshSurfaceOptimizer::smoothSurfaceOptimizer
(
    const labelLongList& selectedPoints
)
{
    this->triMesh();
    updateTriMesh(selectedPoints);

    pointField newPositions(selectedPoints.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 20)
    # endif
    forAll(selectedPoints, i)
    {
        const label bpI = selectedPoints[i];

        newPositions[i] = newPositionSurfaceOptimizer(bpI);
    }

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(newPositions, i)
    {
        const label bpI = selectedPoints[i];

        surfaceModifier.moveBoundaryVertexNoUpdate(bpI, newPositions[i]);
    }

    //- update geometry addressing for moved points
    surfaceModifier.updateGeometry(selectedPoints);
}

bool meshSurfaceOptimizer::untangleSurface
(
    const labelLongList& selectedBoundaryPoints,
    const label nAdditionalLayers
)
{
    Info << "Starting untangling the surface of the volume mesh" << endl;

    bool changed(false);

    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.points();
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();
    surfaceEngine_.boundaryPointEdges();

    if( Pstream::parRun() )
    {
        surfaceEngine_.bpAtProcs();
        surfaceEngine_.globalToLocalBndPointAddressing();
        surfaceEngine_.globalBoundaryPointLabel();
        surfaceEngine_.bpNeiProcs();
    }

    boolList smoothVertex(bPoints.size(), false);
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(selectedBoundaryPoints, i)
        smoothVertex[selectedBoundaryPoints[i]] = true;

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);
    meshSurfaceMapper mapper(surfaceEngine_, meshOctree_);

    bool remapVertex(true);
    label nInvertedTria;
    label nGlobalIter(0);

    labelLongList procBndPoints, movedPoints;
    labelLongList procEdgePoints, movedEdgePoints;

    label minNumInverted(bPoints.size());
    FIFOStack<label> nInvertedHistory;
    pointField minInvertedPoints(bPoints.size());

    do
    {
        label nIter(0), nAfterRefresh(0);

        do
        {
            nInvertedTria =
                findInvertedVertices(smoothVertex, nAdditionalLayers);

            if( nInvertedTria == 0 ) break;

            //- find the min number of inverted points and
            //- add the last number to the stack
            if( nInvertedTria < minNumInverted )
            {
                minNumInverted = nInvertedTria;
                nAfterRefresh = 0;

                # ifdef USE_OMP
                # pragma omp parallel for schedule(dynamic, 100)
                # endif
                forAll(bPoints, bpI)
                    minInvertedPoints[bpI] = points[bPoints[bpI]];
            }

            //- count the number of iteration after the last minimum occurence
            ++nAfterRefresh;

            //- update the stack
            nInvertedHistory.push(nInvertedTria);
            if( nInvertedHistory.size() > 2 )
                nInvertedHistory.pop();

            //- check if the number of inverted points reduces
            bool minimumInStack(false);
            forAllConstIter(FIFOStack<label>, nInvertedHistory, it)
                if( it() == minNumInverted )
                    minimumInStack = true;

            //- stop if the procedure does not minimise
            //- the number of inverted points
            if( !minimumInStack || (nAfterRefresh > 2) )
                break;

            //- find points which will be handled by the smoothers
            changed = true;

            procBndPoints.clear();
            movedPoints.clear();
            procEdgePoints.clear();
            movedEdgePoints.clear();

            forAll(bPoints, bpI)
            {
                if( !smoothVertex[bpI] )
                    continue;

                if( vertexType_[bpI] & PARTITION )
                {
                    movedPoints.append(bpI);

                    if( vertexType_[bpI] & PROCBND )
                        procBndPoints.append(bpI);
                }
                else if( vertexType_[bpI] & EDGE )
                {
                    movedEdgePoints.append(bpI);

                    if( vertexType_[bpI] & PROCBND )
                        procEdgePoints.append(bpI);
                }
            }

            //- smooth edge vertices
            smoothEdgePoints(movedEdgePoints, procEdgePoints);
            if( remapVertex )
                mapper.mapEdgeNodes(movedEdgePoints);
            surfaceModifier.updateGeometry(movedEdgePoints);

            //- use laplacian smoothing
            smoothLaplacianFC(movedPoints, procBndPoints);
            surfaceModifier.updateGeometry(movedPoints);

            //- use surface optimizer
            smoothSurfaceOptimizer(movedPoints);

            if( remapVertex )
                mapper.mapVerticesOntoSurface(movedPoints);

            //- update normals and other geometric data
            surfaceModifier.updateGeometry(movedPoints);

        } while( nInvertedTria && (++nIter < 20) );

        if( nInvertedTria > 0 )
        {
            //- use the combination with the minimu number of inverted points
            meshSurfaceEngineModifier sMod(surfaceEngine_);
            forAll(minInvertedPoints, bpI)
                sMod.moveBoundaryVertexNoUpdate(bpI, minInvertedPoints[bpI]);

            sMod.updateGeometry();
        }

        if( nInvertedTria )
        {
            Info << "Smoothing remaining inverted vertices " << endl;

            movedPoints.clear();
            procBndPoints.clear();
            forAll(smoothVertex, bpI)
                if( smoothVertex[bpI] )
                {
                    movedPoints.append(bpI);

                    if( vertexType_[bpI] & PROCBND )
                        procBndPoints.append(bpI);
                }

            smoothLaplacianFC(movedPoints, procBndPoints, false);

            if( remapVertex )
                mapper.mapVerticesOntoSurface(movedPoints);

            //- update normals and other geometric data
            surfaceModifier.updateGeometry(movedPoints);

            if( nGlobalIter > 3 )
                remapVertex = false;
        }

    } while( nInvertedTria && (++nGlobalIter < 10) );

    Info << "Finished untangling the surface of the volume mesh" << endl;

    return changed;
}

bool meshSurfaceOptimizer::untangleSurface(const label nAdditionalLayers)
{
    labelLongList selectedPts(surfaceEngine_.boundaryPoints().size());
    forAll(selectedPts, i)
        selectedPts[i] = i;

    return untangleSurface(selectedPts, nAdditionalLayers);
}

void meshSurfaceOptimizer::optimizeSurface(const label nIterations)
{
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    //- needed for parallel execution
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryPointEdges();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();
    surfaceEngine_.boundaryPointEdges();

    labelLongList procBndPoints, edgePoints;
    forAll(bPoints, bpI)
    {
        if( vertexType_[bpI] & EDGE )
        {
            edgePoints.append(bpI);

            if( vertexType_[bpI] & PROCBND )
                procBndPoints.append(bpI);
        }
    }

    meshSurfaceMapper mapper(*partitionerPtr_, meshOctree_);

    //- optimize edge vertices
    Info << "Optimizing edges. Iteration:" << flush;
    for(label i=0;i<nIterations;++i)
    {
        Info << "." << flush;

        meshSurfaceEngineModifier bMod(surfaceEngine_);

        smoothEdgePoints(edgePoints, procBndPoints);

        //- project vertices back onto the boundary
        mapper.mapEdgeNodes(edgePoints);

        //- update the geometry information
        bMod.updateGeometry(edgePoints);
    }
    Info << endl;

    //- optimize positions of surface vertices which are not on surface edges
    Info << "Optimizing surface vertices. Iteration:";
    for(label i=0;i<nIterations;++i)
    {
        procBndPoints.clear();

        pointField newPositions(bPoints.size());

        Info << "." << flush;

        meshSurfaceEngineModifier bMod(surfaceEngine_);
        # ifdef USE_OMP
        # pragma omp parallel if( vertexType_.size() > 100 )
        # endif
        {
            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 10)
            # endif
            forAll(bPoints, bpI)
            {
                if( vertexType_[bpI] & PARTITION )
                {
                    if( vertexType_[bpI] & PROCBND )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical
                        # endif
                        {
                            procBndPoints.append(bpI);
                        }

                        continue;
                    }

                    newPositions[bpI] = newPositionLaplacianFC(bpI, true);
                }
            }
        }

        if( Pstream::parRun() )
        {
            nodeDisplacementLaplacianFCParallel(procBndPoints, true);
        }

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 100)
        # endif
        forAll(newPositions, bpI)
        {
            if( vertexType_[bpI] == PARTITION )
                bMod.moveBoundaryVertexNoUpdate(bpI, newPositions[bpI]);
        }

        //- update fields calculated from points
        bMod.updateGeometry();
    }

    Info << endl;

    untangleSurface(0);
}

void meshSurfaceOptimizer::optimizeSurface2D(const label nIterations)
{
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const edgeList& edges = surfaceEngine_.edges();
    const labelList& bp = surfaceEngine_.bp();

    polyMeshGen2DEngine mesh2DEngine
    (
        const_cast<polyMeshGen&>(surfaceEngine_.mesh())
    );
    const boolList& zMinPoint = mesh2DEngine.zMinPoints();

    //- needed for parallel execution
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryPointEdges();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();

    labelLongList procBndPoints, movedPoints, activeEdges, updatePoints;
    forAll(edges, beI)
    {
        const edge& e = edges[beI];

        if( zMinPoint[e.start()] ^ zMinPoint[e.end()] )
        {
            label bpI = bp[e.start()];
            if( !zMinPoint[e.start()] )
                bpI = bp[e.end()];

            if( vertexType_[bpI] & EDGE )
            {
                activeEdges.append(beI);

                updatePoints.append(bp[e.start()]);
                updatePoints.append(bp[e.end()]);

                movedPoints.append(bpI);

                if( vertexType_[bpI] & PROCBND )
                    procBndPoints.append(bpI);
            }
        }
    }

    meshSurfaceMapper2D mapper(surfaceEngine_, meshOctree_);

    //- optimize edge vertices
    meshSurfaceEngineModifier bMod(surfaceEngine_);

    Info << "Optimizing edges. Iteration:" << flush;
    for(label i=0;i<nIterations;++i)
    {
        Info << "." << flush;

        smoothEdgePoints(movedPoints, procBndPoints);

        //- move points with maximum z coordinate
        mesh2DEngine.correctPoints();

        //- map boundary edges to the surface
        mapper.mapVerticesOntoSurfacePatches(activeEdges);

        //- update normal, centres, etc, after the surface has been modified
        bMod.updateGeometry(updatePoints);
    }
    Info << endl;

    //- optimize Pts of surface vertices which are not on surface edges
    procBndPoints.clear();
    movedPoints.clear();
    forAll(bPoints, bpI)
        if( zMinPoint[bPoints[bpI]] && (vertexType_[bpI] & PARTITION) )
        {
            movedPoints.append(bpI);

            if( vertexType_[bpI] & PROCBND )
                procBndPoints.append(bpI);
        }
    Info << "Optimizing surface vertices. Iteration:";
    for(label i=0;i<nIterations;++i)
    {
        Info << "." << flush;

        smoothLaplacianFC(movedPoints, procBndPoints, false);

        //- move the points which are not at minimum z coordinate
        mesh2DEngine.correctPoints();

        //- update geometrical data due to movement of vertices
        bMod.updateGeometry();
    }

    Info << endl;
}

void meshSurfaceOptimizer::untangleSurface2D()
{
    const polyMeshGen& mesh = surfaceEngine_.mesh();
    const faceListPMG& faces = mesh.faces();
    const VRWGraph& pointFaces = mesh.addressingData().pointFaces();

    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const labelList& bp = surfaceEngine_.bp();

    polyMeshGen2DEngine mesh2DEngine(const_cast<polyMeshGen&>(mesh));
    const boolList& zMinPoint = mesh2DEngine.zMinPoints();
    const boolList& activeFace = mesh2DEngine.activeFace();

    //- needed for parallel execution
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryPointEdges();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();

    boolList activeBoundaryPoint(bPoints.size());
    boolList changedFace(activeFace.size(), true);

    label iterationI(0);
    do
    {
        labelHashSet badFaces;
        const label nBadFaces = findBadFaces(badFaces, changedFace);

        Info << "Iteration " << iterationI
             << ". Number of bad faces " << nBadFaces << endl;

        if( nBadFaces == 0 )
            break;

        //- update active points and faces affected by the movement
        //- of active points
        activeBoundaryPoint = false;
        changedFace = false;
        forAllConstIter(labelHashSet, badFaces, it)
        {
            const face& f = faces[it.key()];

            forAll(f, pI)
            {
                if( zMinPoint[f[pI]] )
                {
                    activeBoundaryPoint[bp[f[pI]]] = true;

                    forAllRow(pointFaces, f[pI], pfI)
                        changedFace[pointFaces(f[pI], pfI)] = true;
                }
            }
        }

        if( Pstream::parRun() )
        {
            const Map<label>& globalToLocal =
                surfaceEngine_.globalToLocalBndPointAddressing();
            const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();
            const VRWGraph& bpNeiProcs = surfaceEngine_.bpAtProcs();

            std::map<label, labelLongList> exchangeData;
            forAll(neiProcs, i)
                exchangeData[neiProcs[i]].clear();

            //- collect active points at inter-processor boundaries
            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                if( activeBoundaryPoint[bpI] )
                {
                    forAllRow(bpNeiProcs, bpI, i)
                    {
                        const label neiProc = bpNeiProcs(bpI, i);

                        if( neiProc == Pstream::myProcNo() )
                            continue;

                        exchangeData[neiProc].append(it.key());
                    }
                }
            }

            //- exchange active points among the processors
            labelLongList receivedData;
            help::exchangeMap(exchangeData, receivedData);

            //- ensure that all processors have the same Pts active
            forAll(receivedData, i)
            {
                const label bpI = globalToLocal[receivedData[i]];

                //- activate this boundary point
                activeBoundaryPoint[bpI] = true;

                //- set the changeFaces for the faces attached to this point
                forAllRow(pointFaces, bPoints[bpI], pfI)
                    changedFace[pointFaces(bPoints[bpI], pfI)] = true;
            }
        }

        //- apply smoothing to the activated points
        meshSurfaceEngineModifier bMod(surfaceEngine_);

        labelLongList movedPts, procBndPts, edgePts, procEdgePts;
        forAll(bPoints, bpI)
        {
            if( !activeBoundaryPoint[bpI] )
                continue;

            if( vertexType_[bpI] & EDGE )
            {
                edgePts.append(bpI);

                if( vertexType_[bpI] & PROCBND )
                    procEdgePts.append(bpI);
            }
            else if( vertexType_[bpI] & PARTITION )
            {
                movedPts.append(bpI);

                if( vertexType_[bpI] & PROCBND )
                    procBndPts.append(bpI);
            }
        }

        for(label i=0;i<5;++i)
        {
            smoothEdgePoints(edgePts, procEdgePts);

            bMod.updateGeometry(edgePts);

            smoothSurfaceOptimizer(movedPts);

            bMod.updateGeometry(movedPts);
        }

        //- move the points which are not at minimum z coordinate
        mesh2DEngine.correctPoints();

        //- update geometrical data due to movement of vertices
        bMod.updateGeometry();

        //- update cell centres and face centres
        const_cast<polyMeshGenAddressing&>
        (
            mesh.addressingData()
        ).updateGeometry(changedFace);

    } while( ++iterationI < 20 );

    //- delete invalid data
    mesh.clearAddressingData();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
