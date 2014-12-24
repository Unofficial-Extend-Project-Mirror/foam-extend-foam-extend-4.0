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
#include "plane.H"
#include "meshSurfaceMapper.H"
#include "surfaceOptimizer.H"
#include "refLabelledPoint.H"
#include "helperFunctions.H"

//#define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceOptimizer::nodeDisplacementLaplacianParallel
(
    const labelLongList& nodesToSmooth,
    const bool transformIntoPlane
)
{
    if( returnReduce(nodesToSmooth.size(), sumOp<label>()) == 0 )
        return;

    //- create storage for data
    std::map<label, labelledPoint> localData;

    //- exchange data with other processors
    std::map<label, LongList<refLabelledPoint> > exchangeData;

    const pointField& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const VRWGraph& pPoints = surfaceEngine_.pointPoints();
    const vectorField& pNormals = surfaceEngine_.pointNormals();
    const labelList& globalPointLabel =
        surfaceEngine_.globalBoundaryPointLabel();
    const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
    const Map<label>& globalToLocal =
        surfaceEngine_.globalToLocalBndPointAddressing();

    //- perform smoothing
    forAll(nodesToSmooth, pI)
    {
        const label bpI = nodesToSmooth[pI];

        if( magSqr(pNormals[bpI]) < VSMALL )
            continue;

        const plane pl(points[bPoints[bpI]], pNormals[bpI]);

        //- project points onto the plane
        localData.insert(std::make_pair(bpI, labelledPoint(0, vector::zero)));
        labelledPoint& lpd = localData[bpI];

        forAllRow(pPoints, bpI, ppI)
        {
            const label nei = pPoints(bpI, ppI);

            if( bpAtProcs.sizeOfRow(nei) != 0 )
            {
                label pMin(Pstream::nProcs());
                forAllRow(bpAtProcs, nei, procI)
                {
                    const label procJ = bpAtProcs(nei, procI);
                    if( (procJ < pMin) && bpAtProcs.contains(bpI, procJ) )
                        pMin = procJ;
                }

                if( pMin != Pstream::myProcNo() )
                    continue;
            }

            const point& p = points[bPoints[nei]];
            ++lpd.pointLabel();
            lpd.coordinates() += transformIntoPlane?pl.nearestPoint(p):p;
        }

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            if( exchangeData.find(neiProc) == exchangeData.end() )
            {
                exchangeData.insert
                (
                    std::make_pair(neiProc, LongList<refLabelledPoint>())
                );
            }

            //- add data to the list which will be sent to other processor
            LongList<refLabelledPoint>& dts = exchangeData[neiProc];
            dts.append(refLabelledPoint(globalPointLabel[bpI], lpd));
        }
    }

    //- exchange data with other processors
    LongList<refLabelledPoint> receivedData;
    help::exchangeMap(exchangeData, receivedData);

    forAll(receivedData, i)
    {
        const refLabelledPoint& lp = receivedData[i];
        const label bpI = globalToLocal[lp.objectLabel()];

        labelledPoint& lpd = localData[bpI];

        lpd.pointLabel() += lp.lPoint().pointLabel();
        lpd.coordinates() += lp.lPoint().coordinates();
    }

    forAll(nodesToSmooth, pI)
    {
        const label bpI = nodesToSmooth[pI];

        if( localData.find(bpI) == localData.end() )
            continue;

        //- create new point position
        const labelledPoint& lp = localData[bpI];
        const point newP = lp.coordinates() / lp.pointLabel();

        meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);
        surfaceModifier.moveBoundaryVertexNoUpdate
        (
            bpI,
            newP
        );
    }
}

void meshSurfaceOptimizer::nodeDisplacementLaplacianFCParallel
(
    const labelLongList& nodesToSmooth,
    const bool transformIntoPlane
)
{
    if( returnReduce(nodesToSmooth.size(), sumOp<label>()) == 0 )
        return;

    //- create storage for data
    std::map<label, labelledPoint> localData;

    //- exchange data with other processors
    std::map<label, LongList<refLabelledPoint> > exchangeData;

    const pointField& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const VRWGraph& pFaces = surfaceEngine_.pointFaces();
    const vectorField& faceCentres = surfaceEngine_.faceCentres();
    const vectorField& pNormals = surfaceEngine_.pointNormals();

    const labelList& globalPointLabel =
        surfaceEngine_.globalBoundaryPointLabel();
    const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
    const Map<label>& globalToLocal =
        surfaceEngine_.globalToLocalBndPointAddressing();

    //- perform smoothing
    forAll(nodesToSmooth, pI)
    {
        const label bpI = nodesToSmooth[pI];

        if( magSqr(pNormals[bpI]) < VSMALL )
            continue;

        const plane pl(points[bPoints[bpI]], pNormals[bpI]);

        //- project points onto the plane
        localData.insert(std::make_pair(bpI, labelledPoint(0, vector::zero)));
        labelledPoint& lpd = localData[bpI];

        forAllRow(pFaces, bpI, pfI)
        {
            const point& p = faceCentres[pFaces(bpI, pfI)];
            ++lpd.pointLabel();
            lpd.coordinates() += transformIntoPlane?pl.nearestPoint(p):p;
        }

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            if( exchangeData.find(neiProc) == exchangeData.end() )
            {
                exchangeData.insert
                (
                    std::make_pair(neiProc, LongList<refLabelledPoint>())
                );
            }

            //- add data to the list which will be sent to other processor
            LongList<refLabelledPoint>& dts = exchangeData[neiProc];
            dts.append(refLabelledPoint(globalPointLabel[bpI], lpd));
        }
    }

    //- exchange data with other processors
    LongList<refLabelledPoint> receivedData;
    help::exchangeMap(exchangeData, receivedData);

    forAll(receivedData, i)
    {
        const refLabelledPoint& lp = receivedData[i];
        const label bpI = globalToLocal[lp.objectLabel()];

        labelledPoint& lpd = localData[bpI];

        lpd.pointLabel() += lp.lPoint().pointLabel();
        lpd.coordinates() += lp.lPoint().coordinates();
    }

    pointField newPositions(nodesToSmooth.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 20)
    # endif
    forAll(nodesToSmooth, pI)
    {
        const label bpI = nodesToSmooth[pI];

        if( localData.find(bpI) == localData.end() )
        {
            newPositions[pI] = points[bPoints[bpI]];
            continue;
        }

        //- create new point position
        const labelledPoint& lp = localData[bpI];
        const point newP = lp.coordinates() / lp.pointLabel();

        newPositions[pI] = newP;
    }

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 20)
    # endif
    forAll(newPositions, pI)
    {
        surfaceModifier.moveBoundaryVertexNoUpdate
        (
            nodesToSmooth[pI],
            newPositions[pI]
        );
    }
}

void meshSurfaceOptimizer::nodeDisplacementSurfaceOptimizerParallel
(
    const labelLongList& nodesToSmooth
)
{
    if( returnReduce(nodesToSmooth.size(), sumOp<label>()) == 0 )
        return;

    meshSurfaceEngineModifier surfaceModifier(surfaceEngine_);

    //- exchange data with other processors
    std::map<label, DynList<parTriFace> > m;
    exchangeData(nodesToSmooth, m);

    const pointField& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const vectorField& pNormals = surfaceEngine_.pointNormals();

    //- perform smoothing
    pointField newPositions(nodesToSmooth.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 10)
    # endif
    forAll(nodesToSmooth, pI)
    {
        const label bpI = nodesToSmooth[pI];

        if( magSqr(pNormals[bpI]) < VSMALL )
        {
            newPositions[pI] = points[bPoints[bpI]];
            continue;
        }

        plane pl(points[bPoints[bpI]], pNormals[bpI]);
        point vecX, vecY;
        DynList<point> pts;
        DynList<triFace> trias;

        if( !transformIntoPlaneParallel(bpI, pl, m, vecX, vecY, pts, trias) )
        {
            newPositions[pI] = points[bPoints[bpI]];
            continue;
        }

        surfaceOptimizer so(pts, trias);
        point newPoint = so.optimizePoint(0.001);

        const point newP
        (
            points[bPoints[bpI]] +
            vecX * newPoint.x() +
            vecY * newPoint.y()
        );

        if( help::isnan(newP) || help::isinf(newP) )
        {
            newPositions[pI] = points[bPoints[bpI]];
        }
        else
        {
            newPositions[pI] = newP;
        }
    }

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(newPositions, pI)
        surfaceModifier.moveBoundaryVertexNoUpdate
        (
            nodesToSmooth[pI],
            newPositions[pI]
        );

    //- make sure that moved points have the same coordinates on all processors
    surfaceModifier.syncVerticesAtParallelBoundaries(nodesToSmooth);
}

void meshSurfaceOptimizer::edgeNodeDisplacementParallel
(
    const labelLongList& nodesToSmooth
)
{
    std::map<label, DynList<labelledPoint, 2> > mPts;

    //- insert labels into the map
    const labelList& globalPointLabel =
        surfaceEngine_.globalBoundaryPointLabel();

    forAll(nodesToSmooth, nI)
    {
        mPts.insert
        (
            std::make_pair
            (
                globalPointLabel[nodesToSmooth[nI]],
                DynList<labelledPoint, 2>()
            )
        );
    }

    //- add local edge points
    const pointField& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const edgeList& edges = surfaceEngine_.edges();
    const VRWGraph& bpEdges = surfaceEngine_.boundaryPointEdges();
    const  labelList& bp = surfaceEngine_.bp();

    forAll(nodesToSmooth, nI)
    {
        const label bpI = nodesToSmooth[nI];
        DynList<labelledPoint, 2>& neiPoints = mPts[globalPointLabel[bpI]];

        forAllRow(bpEdges, bpI, epI)
        {
            const edge& e = edges[bpEdges(bpI, epI)];
            const label pI = bp[e.otherVertex(bPoints[bpI])];
            if( vertexType_[pI] & (EDGE+CORNER) )
            {
                neiPoints.append
                (
                    labelledPoint(globalPointLabel[pI], points[bPoints[pI]])
                );
            }
        }
    }

    //- start preparing data which can be exchanged with other processors
    const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();
    const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();

    std::map<label, LongList<refLabelledPoint> > mProcs;
    forAll(neiProcs, procI)
    {
        const label neiProc = neiProcs[procI];

        if( neiProc == Pstream::myProcNo() )
            continue;

        mProcs.insert(std::make_pair(neiProc, LongList<refLabelledPoint>()));
    }

    forAll(nodesToSmooth, nI)
    {
        const label bpI = nodesToSmooth[nI];

        const DynList<labelledPoint, 2>& neiPoints =
            mPts[globalPointLabel[bpI]];

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);

            if( neiProc == Pstream::myProcNo() )
                continue;

            forAll(neiPoints, npI)
            {
                LongList<refLabelledPoint>& neiProcPts = mProcs[neiProc];
                neiProcPts.append
                (
                    refLabelledPoint(globalPointLabel[bpI], neiPoints[npI])
                );
            }
        }
    }

    //- exchange data with other processors
    LongList<refLabelledPoint> receivedData;
    help::exchangeMap(mProcs, receivedData);

    forAll(receivedData, prI)
    {
        const refLabelledPoint& lp = receivedData[prI];
        DynList<labelledPoint, 2>& lPts = mPts[lp.objectLabel()];
        lPts.appendIfNotIn(receivedData[prI].lPoint());
    }

    //- Finally, the data is ready to start smoothing
    meshSurfaceEngineModifier sm(surfaceEngine_);

    pointField newPositions(nodesToSmooth.size());
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 20)
    # endif
    forAll(nodesToSmooth, pI)
    {
        const label bpI = nodesToSmooth[pI];

        const DynList<labelledPoint, 2>& nPts = mPts[globalPointLabel[bpI]];
        point newP(vector::zero);
        forAll(nPts, ppI)
            newP += nPts[ppI].coordinates();

        if( nPts.size() == 2 )
        {
            newP /= nPts.size();
            newPositions[pI] = newP;
        }
        else
        {
            newPositions[pI] = points[bPoints[bpI]];
        }
    }

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 20)
    # endif
    forAll(newPositions, pI)
        sm.moveBoundaryVertexNoUpdate(nodesToSmooth[pI], newPositions[pI]);
}

void meshSurfaceOptimizer::exchangeData
(
    const labelLongList& nodesToSmooth,
    std::map<label, DynList<parTriFace> >& m
) const
{
    /*
    const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
    const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();

    std::map<label, LongList<parTriFace> > shareData;
    forAll(neiProcs, procI)
    {
        const label neiProc = neiProcs[procI];

        if( neiProc == Pstream::myProcNo() )
            continue;

        shareData.insert(std::make_pair(neiProc, LongList<parTriFace>()));
    }

    //- create data which will be sent to other processors
    const pointField& points = surfaceEngine_.points();
    const labelList& globalPointLabel =
        surfaceEngine_.globalBoundaryPointLabel();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const LongList<triFace>& triangles = this->triangles();
    const VRWGraph& pTriangles = pointTriangles();

    m.clear();
    std::map<label, DynList<parTriFace> >::iterator pIter;
    forAll(nodesToSmooth, nI)
    {
        const label bpI = nodesToSmooth[nI];

        pIter = m.find(globalPointLabel[bpI]);
        if( pIter == m.end() )
        {
            m.insert
            (
                std::make_pair(globalPointLabel[bpI], DynList<parTriFace>())
            );

            pIter = m.find(globalPointLabel[bpI]);
        }

        forAll(pTriangles[bpI], ptI)
        {
            const triFace& tri = triangles[pTriangles[bpI][ptI]];

            parTriFace ptf
            (
                globalPointLabel[tri[0]],
                globalPointLabel[tri[1]],
                globalPointLabel[tri[2]],
                triangle<point, point>
                (
                    points[bPoints[tri[0]]],
                    points[bPoints[tri[1]]],
                    points[bPoints[tri[2]]]
                )
            );

            pIter->second.append(ptf);
        }

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            LongList<parTriFace>& shData = shareData[neiProc];

            const DynList<parTriFace>& localData = pIter->second;

            forAll(localData, i)
                shData.append(localData[i]);
        }
    }

    //- send data to other processors
    LongList<parTriFace> receivedData;
    help::exchangeMap(shareData, receivedData);

    forAll(receivedData, tI)
    {
        const label gpI = receivedData[tI].globalLabelOfPoint(0);

        pIter = m.find(gpI);
        if( pIter == m.end() )
        {
            FatalErrorIn
            (
                "void meshSurfaceOptimizer::exchangeData("
                "const labelLongList& nodesToSmooth,"
                "map<label, DynList<parTriFace> >& m"
                ") const"
            ) << "Unknown point " << gpI << abort(FatalError);
        }

        pIter->second.append(receivedData[tI]);
    }
    */
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
