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

#include "meshOctree.H"
#include "triSurf.H"
#include "triSurfacePartitioner.H"
#include "meshSurfaceMapper2D.H"
#include "meshSurfaceEngine.H"
#include "polyMeshGen2DEngine.H"
#include "meshSurfaceEngineModifier.H"
#include "meshSurfacePartitioner.H"
#include "labelledScalar.H"
#include "polyMeshGenAddressing.H"
#include "boundBox.H"

#include "helperFunctionsPar.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGMapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceMapper2D::findMappingDistance
(
    const labelLongList& edgesToMap,
    std::map<label, scalar>& mappingDistance
) const
{
    const vectorField& faceCentres = surfaceEngine_.faceCentres();
    const VRWGraph& eFaces = surfaceEngine_.edgeFaces();
    const edgeList& edges = surfaceEngine_.edges();
    const pointFieldPMG& points = surfaceEngine_.points();

    //- generate search distance for corner edges
    mappingDistance.clear();
    forAll(edgesToMap, i)
    {
        const label beI = edgesToMap[i];

        mappingDistance.insert(std::make_pair(beI, 0.0));
        std::map<label, scalar>::iterator mIter = mappingDistance.find(beI);

        const point p = edges[beI].centre(points);
        forAllRow(eFaces, beI, efI)
        {
            const scalar d = magSqr(faceCentres[eFaces(beI, efI)] - p);
            mIter->second = Foam::max(mIter->second, d);
        }

        //- safety factor
        mIter->second *= 16.0;
    }

    if( Pstream::parRun() )
    {
        //- make sure that corner edges at parallel boundaries
        //- have the same range in which they accept the corners
        const VRWGraph& beAtProcs = surfaceEngine_.beAtProcs();
        const labelList& globalEdgeLabel =
            surfaceEngine_.globalBoundaryEdgeLabel();
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndEdgeAddressing();

        //- create the map for exchanging data
        std::map<label, DynList<labelledScalar> > exchangeData;
        const DynList<label>& neiProcs = surfaceEngine_.beNeiProcs();
        forAll(neiProcs, i)
            exchangeData.insert
            (
                std::make_pair(neiProcs[i], DynList<labelledScalar>()));

        forAll(edgesToMap, eI)
        {
            const label beI = edgesToMap[eI];

            forAllRow(beAtProcs, beI, i)
            {
                const label neiProc = beAtProcs(beI, i);
                if( neiProc == Pstream::myProcNo() )
                    continue;

                exchangeData[neiProc].append
                (
                    labelledScalar(globalEdgeLabel[beI], mappingDistance[beI])
                );
            }
        }

        //- exchange data between processors
        LongList<labelledScalar> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        //- select the maximum mapping distance for processor points
        forAll(receivedData, i)
        {
            const labelledScalar& ls = receivedData[i];

            const label beI = globalToLocal[ls.scalarLabel()];

            //- choose the maximum value for the mapping distance
            std::map<label, scalar>::iterator mIter = mappingDistance.find(beI);
            mIter->second = Foam::max(mIter->second, ls.value());
        }
    }
}

void meshSurfaceMapper2D::mapToSmallestDistance(LongList<parMapperHelper>& parE)
{
    if( !Pstream::parRun() )
        return;

    std::map<label, LongList<parMapperHelper> > exchangeData;
    const DynList<label>& neiProcs = surfaceEngine_.beNeiProcs();
    forAll(neiProcs, i)
        exchangeData.insert
        (
            std::make_pair(neiProcs[i], LongList<parMapperHelper>())
        );

    const labelList& bp = surfaceEngine_.bp();
    const edgeList& edges = surfaceEngine_.edges();
    const VRWGraph& beAtProcs = surfaceEngine_.beAtProcs();
    const labelList& globalEdgeLabel =
        surfaceEngine_.globalBoundaryEdgeLabel();
    const Map<label>& globalToLocal =
        surfaceEngine_.globalToLocalBndEdgeAddressing();
    const pointFieldPMG& points = surfaceEngine_.points();

    Map<label> beToList(parE.size());

    forAll(parE, i)
    {
        const label beI = parE[i].globalLabel();
        beToList.insert(beI, i);

        forAllRow(beAtProcs, beI, procI)
        {
            const label neiProc = beAtProcs(beI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            exchangeData[neiProc].append
            (
                parMapperHelper
                (
                    parE[i].coordinates(),
                    parE[i].movingDistance(),
                    globalEdgeLabel[beI],
                    parE[i].pointPatch()
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

        const label beI = globalToLocal[ph.globalLabel()];

        parMapperHelper& phOrig = parE[beToList[beI]];
        if( phOrig.movingDistance() < ph.movingDistance() )
        {
            const edge& e = edges[beI];

            point newP = ph.coordinates();

            newP.z() = points[e.start()].z();
            surfModifier.moveBoundaryVertex(bp[e.start()], newP);

            newP.z() = points[e.end()].z();
            surfModifier.moveBoundaryVertex(bp[e.end()], newP);

            phOrig = ph;
        }
    }

    surfModifier.updateVertexNormals();
}

void meshSurfaceMapper2D::adjustZCoordinates()
{
    //- create the bounding box of the surface mesh
    const boundBox sbb(meshOctree_.surface().points());

    const labelList& bp = surfaceEngine_.bp();
    const pointFieldPMG& points = surfaceEngine_.points();

    meshSurfaceEngineModifier modifier(surfaceEngine_);

    const polyMeshGen2DEngine& mesh2DEngine = this->mesh2DEngine();
    const boolList& minZPoint = mesh2DEngine.zMinPoints();
    const boolList& maxZPoint = mesh2DEngine.zMaxPoints();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(minZPoint, pointI)
    {
        point newP = points[pointI];

        if( minZPoint[pointI] )
        {
            newP.z() = sbb.min().z();
        }
        else if( maxZPoint[pointI] )
        {
            newP.z() = sbb.max().z();
        }
        else
        {
            FatalErrorIn
            (
                "void meshSurfaceMapper2D::adjustZCoordinates()"
            ) << "This mesh is not in the x-y plane!" << exit(FatalError);
        }

        modifier.moveBoundaryVertexNoUpdate(bp[pointI], newP);
    }

    modifier.updateGeometry();
}

void meshSurfaceMapper2D::mapVerticesOntoSurface()
{
    labelLongList edgesToMap;

    forAll(activeBoundaryEdges_, beI)
        edgesToMap.append(activeBoundaryEdges_[beI]);

    mapVerticesOntoSurface(edgesToMap);
}

void meshSurfaceMapper2D::mapVerticesOntoSurface(const labelLongList& edgesToMap)
{
    const edgeList& edges = surfaceEngine_.edges();
    const labelList& bp = surfaceEngine_.bp();
    const pointFieldPMG& points = surfaceEngine_.points();

    const VRWGraph* beAtProcsPtr(NULL);
    if( Pstream::parRun() )
        beAtProcsPtr = &surfaceEngine_.beAtProcs();

    labelLongList nodesToMap;
    forAll(edgesToMap, i)
    {
        const edge& e = edges[edgesToMap[i]];
        nodesToMap.append(bp[e.start()]);
        nodesToMap.append(bp[e.end()]);
    }

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);
    LongList<parMapperHelper> parallelBndEdges;

    # ifdef USE_OMP
    const label size = edgesToMap.size();
    # pragma omp parallel for if( size > 1000 ) shared(parallelBndEdges) \
    schedule(dynamic, Foam::max(1, size / (3 * omp_get_max_threads())))
    # endif
    forAll(edgesToMap, i)
    {
        const label beI = edgesToMap[i];
        const edge& e = edges[beI];

        # ifdef DEBUGMapping
        Info << nl << "Mapping edge " << beI << " with nodes "
             << e << endl;
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
            points[e.start()]
        );

        mapPoint.z() = points[e.start()].z();
        surfaceModifier.moveBoundaryVertexNoUpdate(bp[e.start()], mapPoint);

        mapPoint.z() = points[e.end()].z();
        surfaceModifier.moveBoundaryVertexNoUpdate(bp[e.end()], mapPoint);

        if( beAtProcsPtr && beAtProcsPtr->sizeOfRow(beI) )
        {
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            parallelBndEdges.append
            (
                parMapperHelper
                (
                    mapPoint,
                    dSq,
                    beI,
                    patch
                )
            );
        }

        # ifdef DEBUGMapping
        Info << "Mapped edge " << edges[beI] << endl;
        # endif
    }

    surfaceModifier.updateGeometry(nodesToMap);

    mapToSmallestDistance(parallelBndEdges);
}

void meshSurfaceMapper2D::mapCorners()
{
    const edgeList& edges = surfaceEngine_.edges();

    const meshSurfacePartitioner& mPart = meshPartitioner();
    const labelHashSet& cornerPoints = mPart.corners();

    labelLongList selectedEdges;
    forAll(activeBoundaryEdges_, eI)
    {
        const edge& e = edges[activeBoundaryEdges_[eI]];

        if( cornerPoints.found(e.start()) || cornerPoints.found(e.end()) )
            selectedEdges.append(activeBoundaryEdges_[eI]);
    }

    mapCorners(selectedEdges);
}

void meshSurfaceMapper2D::mapCorners(const labelLongList& edgesToMap)
{
    const meshSurfacePartitioner& mPart = meshPartitioner();
    const labelHashSet& corners = mPart.corners();
    const VRWGraph& pPatches = mPart.pointPatches();

    const pointFieldPMG& points = surfaceEngine_.points();
    const edgeList& edges = surfaceEngine_.edges();
    const labelList& bp = surfaceEngine_.bp();

    std::map<label, scalar> mappingDistance;
    findMappingDistance(edgesToMap, mappingDistance);

    //- for every corner in the mesh surface find the nearest corner in the
    //- triSurface
    meshSurfaceEngineModifier sMod(surfaceEngine_);

    # ifdef USE_OMP
    # pragma omp parallel for if( edgesToMap.size() > 10 )
    # endif
    forAll(edgesToMap, eI)
    {
        const label beI = edgesToMap[eI];
        const edge& be = edges[beI];

        const label bps = bp[be.start()];
        const label bpe = bp[be.end()];

        DynList<label> ePatches;
        forAllRow(pPatches, bps, i)
        {
            if( pPatches.contains(bpe, pPatches(bps, i)) )
                ePatches.append(pPatches(bps, i));
        }

        if( !corners.found(bps) || !corners.found(bpe) )
            FatalErrorIn
            (
                "meshSurfaceMapper2D::mapCorners(const labelLongList&)"
            ) << "Trying to map a point that is not a corner"
                << abort(FatalError);

        const point& p = points[be.start()];
        const scalar maxDist = mappingDistance[beI];

        //- find the nearest position to the given point patches
        point mapPointApprox(p);
        scalar distSqApprox;

        label iter(0);
        while( iter++ < 40 )
        {
            point newP(vector::zero);
            forAll(ePatches, epI)
            {
                point np;
                label nt;
                meshOctree_.findNearestSurfacePointInRegion
                (
                    np,
                    distSqApprox,
                    nt,
                    ePatches[epI],
                    mapPointApprox
                );

                newP += np;
            }

            newP /= ePatches.size();
            if( magSqr(newP - mapPointApprox) < 1e-8 * maxDist )
                break;

            mapPointApprox = newP;
        }

        distSqApprox = magSqr(mapPointApprox - p);

        //- find the nearest triSurface corner for the given corner
        scalar distSq;
        point mapPoint;
        label nse;
        meshOctree_.findNearestEdgePoint(mapPoint, distSq, nse, p, ePatches);

        if( distSq > mappingDistance[beI] )
        {
            const vector disp = mapPoint - p;

            mapPoint = Foam::sqrt(mappingDistance[beI] / distSq) * disp + p;
            distSq = mappingDistance[beI];
        }

        if( distSq > 1.2 * distSqApprox )
        {
            mapPoint = mapPointApprox;
        }

        //- move the point to the nearest corner
        mapPoint.z() = points[be.start()].z();
        sMod.moveBoundaryVertexNoUpdate(bps, mapPoint);

        mapPoint.z() = points[be.end()].z();
        sMod.moveBoundaryVertexNoUpdate(bpe, mapPoint);
    }

    labelLongList nodesToMap;
    forAll(edgesToMap, eI)
    {
        const edge& be = edges[edgesToMap[eI]];

        nodesToMap.append(bp[be.start()]);
        nodesToMap.append(bp[be.end()]);
    }

    sMod.updateGeometry(nodesToMap);
}

void meshSurfaceMapper2D::mapVerticesOntoSurfacePatches()
{
    labelLongList edgesToMap;

    forAll(activeBoundaryEdges_, beI)
        edgesToMap.append(activeBoundaryEdges_[beI]);

    mapVerticesOntoSurfacePatches(edgesToMap);
}

void meshSurfaceMapper2D::mapVerticesOntoSurfacePatches
(
    const labelLongList& edgesToMap
)
{
    //- map the selected edges on their patches
    const edgeList& edges = surfaceEngine_.edges();
    const pointFieldPMG& points = surfaceEngine_.points();
    const VRWGraph& edgePatches = surfaceEngine_.edgePatches();
    const labelList& bp = surfaceEngine_.bp();

    const VRWGraph* beAtProcsPtr(NULL);
    if( Pstream::parRun() )
        beAtProcsPtr = &surfaceEngine_.beAtProcs();

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);
    LongList<parMapperHelper> parallelBndEdges;
    labelLongList selectedCorners;

    # ifdef USE_OMP
    const label size = edgesToMap.size();
    # pragma omp parallel for if( size > 1000 ) shared(parallelBndEdges) \
    schedule(dynamic, Foam::max(50, size / (3 * omp_get_max_threads())))
    # endif
    forAll(edgesToMap, nI)
    {
        const label beI = edgesToMap[nI];

        if( edgePatches.sizeOfRow(beI) == 2 )
        {
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            selectedCorners.append(beI);

            continue;
        }

        const edge& e = edges[beI];

        const point& p = points[e.start()];
        point mapPoint;
        scalar dSq;
        label nt;

        meshOctree_.findNearestSurfacePointInRegion
        (
            mapPoint,
            dSq,
            nt,
            edgePatches(beI, 0),
            p
        );

        mapPoint.z() = points[e.start()].z();
        surfaceModifier.moveBoundaryVertexNoUpdate(bp[e.start()], mapPoint);

        mapPoint.z() = points[e.end()].z();
        surfaceModifier.moveBoundaryVertexNoUpdate(bp[e.end()], mapPoint);

        if( beAtProcsPtr && beAtProcsPtr->sizeOfRow(beI) )
        {
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            parallelBndEdges.append
            (
                parMapperHelper
                (
                    mapPoint,
                    dSq,
                    beI,
                    -1
                )
            );
        }

        # ifdef DEBUGMapping
        Info << "Mapped edge " << edges[beI] << endl;
        # endif
    }

    mapToSmallestDistance(parallelBndEdges);

    mapCorners(selectedCorners);

    //- update geometry at moved vertices
    selectedCorners.clear();
    forAll(edgesToMap, eI)
    {
        const edge& e = edges[edgesToMap[eI]];

        selectedCorners.append(bp[e.start()]);
        selectedCorners.append(bp[e.end()]);
    }

    surfaceModifier.updateGeometry(selectedCorners);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
