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
#include "meshSurfaceEngineModifier.H"
#include "meshSurfaceMapper.H"
#include "meshSurfacePartitioner.H"
#include "meshOctree.H"
#include "triSurf.H"
#include "helperFunctionsPar.H"
#include "helperFunctions.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGMapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Private member functions

void meshSurfaceMapper::selectNodesAtParallelBnd(const labelLongList& selNodes)
{
    if( !Pstream::parRun() )
        return;

    std::map<label, labelLongList> exchangeData;
    const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();
    forAll(neiProcs, i)
        exchangeData.insert(std::make_pair(neiProcs[i], labelLongList()));

    const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
    const labelList& globalPointLabel =
        surfaceEngine_.globalBoundaryPointLabel();
    const Map<label>& globalToLocal =
        surfaceEngine_.globalToLocalBndPointAddressing();

    boolList selectedNode(bpAtProcs.size(), false);

    forAll(selNodes, i)
    {
        const label bpI = selNodes[i];

        selectedNode[bpI] = true;

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            exchangeData[neiProc].append(globalPointLabel[bpI]);
        }
    }

    //- exchange data
    labelLongList receivedData;
    help::exchangeMap(exchangeData, receivedData);

    forAll(receivedData, i)
    {
        if( !selectedNode[globalToLocal[receivedData[i]]] )
        {
            selectedNode[globalToLocal[receivedData[i]]] = true;
            const_cast<labelLongList&>(selNodes).append
            (
                globalToLocal[receivedData[i]]
            );
        }
    }
}

void meshSurfaceMapper::mapToSmallestDistance(LongList<parMapperHelper>& parN)
{
    if( !Pstream::parRun() )
        return;

    std::map<label, LongList<parMapperHelper> > exchangeData;
    const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();
    forAll(neiProcs, i)
        exchangeData.insert
        (
            std::make_pair(neiProcs[i], LongList<parMapperHelper>())
        );

    const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
    const labelList& globalPointLabel =
        surfaceEngine_.globalBoundaryPointLabel();
    const Map<label>& globalToLocal =
        surfaceEngine_.globalToLocalBndPointAddressing();

    Map<label> bpToList(parN.size());

    forAll(parN, i)
    {
        const label bpI = parN[i].globalLabel();
        bpToList.insert(bpI, i);

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            exchangeData[neiProc].append
            (
                parMapperHelper
                (
                    parN[i].coordinates(),
                    parN[i].movingDistance(),
                    globalPointLabel[bpI],
                    parN[i].pointPatch()
                )
            );
        }
    }

    //- exchange data
    LongList<parMapperHelper> receivedData;
    help::exchangeMap(exchangeData, receivedData);

    //- select the point with the smallest moving distance
    meshSurfaceEngineModifier surfModifier(surfaceEngine_);
    forAll(receivedData, i)
    {
        const parMapperHelper& ph = receivedData[i];

        const label bpI = globalToLocal[ph.globalLabel()];

        parMapperHelper& phOrig = parN[bpToList[bpI]];
        if( phOrig.movingDistance() < ph.movingDistance() )
        {
            surfModifier.moveBoundaryVertexNoUpdate(bpI, ph.coordinates());
            phOrig = ph;
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceMapper::mapNodeToPatch(const label bpI, const label patchI)
{
    label patch, nt;
    point mapPoint;
    scalar dSq;

    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const point p = points[bPoints[bpI]];

    if( patchI < 0 )
    {
        meshOctree_.findNearestSurfacePoint(mapPoint, dSq, nt, patch, p);
    }
    else
    {
        meshOctree_.findNearestSurfacePointInRegion
        (
            mapPoint,
            dSq,
            nt,
            patchI,
            p
        );
    }

    meshSurfaceEngineModifier surfModifier(surfaceEngine_);
    surfModifier.moveBoundaryVertex(bpI, mapPoint);
}

void meshSurfaceMapper::mapVerticesOntoSurface()
{
    Info << "Mapping vertices onto surface" << endl;

    labelLongList nodesToMap(surfaceEngine_.boundaryPoints().size());
    forAll(nodesToMap, i)
        nodesToMap[i] = i;

    mapVerticesOntoSurface(nodesToMap);

    Info << "Finished mapping vertices onto surface" << endl;
}

void meshSurfaceMapper::mapVerticesOntoSurface(const labelLongList& nodesToMap)
{
    const labelList& boundaryPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.points();

    const VRWGraph* bpAtProcsPtr(NULL);
    if( Pstream::parRun() )
        bpAtProcsPtr = &surfaceEngine_.bpAtProcs();

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);
    LongList<parMapperHelper> parallelBndNodes;

    # ifdef USE_OMP
    const label size = nodesToMap.size();
    # pragma omp parallel for if( size > 1000 ) shared(parallelBndNodes) \
    schedule(dynamic, Foam::max(1, size / (3 * omp_get_max_threads())))
    # endif
    forAll(nodesToMap, i)
    {
        const label bpI = nodesToMap[i];

        # ifdef DEBUGMapping
        Info << nl << "Mapping vertex " << bpI << " with coordinates "
            << points[boundaryPoints[bpI]] << endl;
        # endif

        label patch, nt;
        point mapPoint;
        scalar dSq;

        meshOctree_.findNearestSurfacePoint
        (
            mapPoint,
            dSq,
            nt,
            patch,
            points[boundaryPoints[bpI]]
        );

        surfaceModifier.moveBoundaryVertexNoUpdate(bpI, mapPoint);

        if( bpAtProcsPtr && bpAtProcsPtr->sizeOfRow(bpI) )
        {
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            parallelBndNodes.append
            (
                parMapperHelper
                (
                    mapPoint,
                    dSq,
                    bpI,
                    patch
                )
            );
        }

        # ifdef DEBUGMapping
        Info << "Mapped point " << points[boundaryPoints[bpI]] << endl;
        # endif
    }

    //- make sure that the points are at the nearest location on the surface
    mapToSmallestDistance(parallelBndNodes);

    //- re-calculate face normals, point normals, etc.
    surfaceModifier.updateGeometry(nodesToMap);
}

void meshSurfaceMapper::mapVerticesOntoSurfacePatches()
{
    Info << "Mapping vertices with respect to surface patches" << endl;

    labelLongList nodesToMap(surfaceEngine_.boundaryPoints().size());
    forAll(nodesToMap, i)
        nodesToMap[i] = i;

    mapVerticesOntoSurfacePatches(nodesToMap);

    Info << "Finished mapping vertices with respect to surface patches" << endl;
}

void meshSurfaceMapper::mapVerticesOntoSurfacePatches
(
    const labelLongList& nodesToMap
)
{
    const meshSurfacePartitioner& mPart = meshPartitioner();
    const labelHashSet& cornerPoints = mPart.corners();
    const labelHashSet& edgePoints = mPart.edgePoints();
    const VRWGraph& pointPatches = mPart.pointPatches();

    boolList treatedPoint(surfaceEngine_.boundaryPoints().size(), false);

    //- find corner and edge points
    labelLongList selectedCorners, selectedEdges;
    forAll(nodesToMap, i)
    {
        if( cornerPoints.found(nodesToMap[i]) )
        {
            treatedPoint[nodesToMap[i]] = true;
            selectedCorners.append(nodesToMap[i]);
        }
        else if( edgePoints.found(nodesToMap[i]) )
        {
            treatedPoint[nodesToMap[i]] = true;
            selectedEdges.append(nodesToMap[i]);
        }
    }

    //- map the remaining selected points
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.points();

    const VRWGraph* bpAtProcsPtr(NULL);
    if( Pstream::parRun() )
        bpAtProcsPtr = &surfaceEngine_.bpAtProcs();

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);
    LongList<parMapperHelper> parallelBndNodes;

    # ifdef USE_OMP
    const label size = nodesToMap.size();
    # pragma omp parallel for if( size > 1000 ) shared(parallelBndNodes) \
    schedule(dynamic, Foam::max(1, size / (3 * omp_get_max_threads())))
    # endif
    forAll(nodesToMap, nI)
    {
        const label bpI = nodesToMap[nI];

        if( treatedPoint[bpI] )
            continue;

        const point& p = points[bPoints[bpI]];
        point mapPoint;
        scalar dSq;
        label nt;
        meshOctree_.findNearestSurfacePointInRegion
        (
            mapPoint,
            dSq,
            nt,
            pointPatches(bpI, 0),
            p
        );

        surfaceModifier.moveBoundaryVertexNoUpdate(bpI, mapPoint);

        if( bpAtProcsPtr && bpAtProcsPtr->sizeOfRow(bpI) )
        {
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            {
                parallelBndNodes.append
                (
                    parMapperHelper
                    (
                        mapPoint,
                        dSq,
                        bpI,
                        -1
                    )
                );
            }
        }

        # ifdef DEBUGMapping
        Info << "Mapped point " << points[boundaryPoints[bpI]] << endl;
        # endif
    }

    //- map vertices at inter-processor boundaries to the nearest location
    //- on the surface
    mapToSmallestDistance(parallelBndNodes);

    //- update face normals, point normals, etc.
    surfaceModifier.updateGeometry(nodesToMap);

    //- map edge nodes
    mapEdgeNodes(selectedEdges);

    //- map corner vertices
    mapCorners(selectedCorners);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
