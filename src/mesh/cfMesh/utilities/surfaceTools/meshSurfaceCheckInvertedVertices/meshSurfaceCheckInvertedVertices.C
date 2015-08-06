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

#include "meshSurfaceCheckInvertedVertices.H"
#include "meshSurfacePartitioner.H"
#include "boolList.H"
#include "demandDrivenData.H"
#include "refLabelledPoint.H"
#include "helperFunctions.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceCheckInvertedVertices::checkVertices()
{
    const labelList& facePatch = surfacePartitioner_.boundaryFacePatches();
    const meshSurfaceEngine& mse = surfacePartitioner_.surfaceEngine();
    const pointFieldPMG& points = mse.points();
    const labelList& bp = mse.bp();
    const VRWGraph& pointFaces = mse.pointFaces();
    const VRWGraph& pointInFaces = mse.pointInFaces();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const vectorField& pNormals = mse.pointNormals();
    const vectorField& fCentres = mse.faceCentres();
    const vectorField& fNormals = mse.faceNormals();

    const labelHashSet& corners = surfacePartitioner_.corners();
    const labelHashSet& edgePoints = surfacePartitioner_.edgePoints();

    typedef std::map<label, vector> ltvMap;
    typedef std::map<label, ltvMap> lltvMap;
    lltvMap pointPatchNormal;

    forAllConstIter(labelHashSet, corners, it)
    {
        const label bpI = it.key();

        if( activePointsPtr_ && !activePointsPtr_->operator[](bpI))
            continue;

        ltvMap& patchNormal = pointPatchNormal[bpI];

        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            const label patchI = facePatch[bfI];

            if( patchNormal.find(patchI) == patchNormal.end() )
            {
                patchNormal[patchI] = fNormals[bfI];
            }
            else
            {
                patchNormal[patchI] += fNormals[bfI];
            }
        }
    }

    forAllConstIter(labelHashSet, edgePoints, it)
    {
        const label bpI = it.key();

        if( activePointsPtr_ && !activePointsPtr_->operator[](bpI))
            continue;

        ltvMap& patchNormal = pointPatchNormal[bpI];

        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            const label patchI = facePatch[bfI];

            if( patchNormal.find(patchI) == patchNormal.end() )
            {
                patchNormal[patchI] = fNormals[bfI];
            }
            else
            {
                patchNormal[patchI] += fNormals[bfI];
            }
        }
    }

    if( Pstream::parRun() )
    {
        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();
        const DynList<label>& neiProcs = mse.bpNeiProcs();
        const VRWGraph& bpAtProcs = mse.bpAtProcs();

        std::map<label, LongList<refLabelledPoint> > exchangeData;
        forAll(neiProcs, i)
            exchangeData[neiProcs[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            if( pointPatchNormal.find(bpI) != pointPatchNormal.end() )
            {
                const ltvMap& patchNormal = pointPatchNormal[bpI];

                forAllRow(bpAtProcs, bpI, i)
                {
                    const label neiProc = bpAtProcs(bpI, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    forAllConstIter(ltvMap, patchNormal, pIt)
                        exchangeData[neiProc].append
                        (
                            refLabelledPoint
                            (
                                it.key(),
                                labelledPoint(pIt->first, pIt->second)
                            )
                        );
                }
            }
        }

        LongList<refLabelledPoint> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const refLabelledPoint& rlp = receivedData[i];
            const label bpI = globalToLocal[rlp.objectLabel()];

            ltvMap& patchNormal = pointPatchNormal[bpI];

            const labelledPoint& lp = rlp.lPoint();
            patchNormal[lp.pointLabel()] += lp.coordinates();
        }
    }

    forAllIter(lltvMap, pointPatchNormal, it)
    {
        ltvMap& patchNormal = pointPatchNormal[it->first];

        forAllIter(ltvMap, patchNormal, pIt)
        {
            const scalar magv = mag(pIt->second) + VSMALL;

            pIt->second /= magv;
        }
    }

    invertedVertices_.clear();

    # ifdef USE_OMP
    # pragma omp parallel for if( pointFaces.size() > 100 ) \
    schedule(dynamic, 20)
    # endif
    forAll(pointFaces, bpI)
    {
        if( activePointsPtr_ && !activePointsPtr_->operator[](bpI) )
            continue;

        forAllRow(pointFaces, bpI, pfI)
        {
            const label pI = pointInFaces(bpI, pfI);
            const label bfI = pointFaces(bpI, pfI);

            vector pNormal = pNormals[bpI];

            if( pointPatchNormal.find(bpI) != pointPatchNormal.end() )
                pNormal = pointPatchNormal[bpI][facePatch[bfI]];

            const face& bf = bFaces[bfI];

            //- chech the first triangle (with the next node)
            triangle<point, point> triNext
            (
                points[bf[pI]],
                points[bf.nextLabel(pI)],
                fCentres[bfI]
            );

            vector nNext = triNext.normal();
            scalar mNext = mag(nNext);

            //- face has zero area
            if( mNext < VSMALL )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                invertedVertices_.insert(bf[pI]);

                break;
            }
            else
            {
                nNext /= mNext;
            }

            //- collocated points
            if( magSqr(triNext.a() - triNext.b()) < VSMALL )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                invertedVertices_.insert(bf[pI]);

                break;
            }
            if( magSqr(triNext.c() - triNext.a()) < VSMALL )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                invertedVertices_.insert(bf[pI]);

                break;
            }

            //- normal vector is not visible
            if( (nNext & pNormal) < 0.0 )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                invertedVertices_.insert(bf[pI]);

                break;
            }

            //- check the second triangle (with previous node)
            triangle<point, point> triPrev
            (
                points[bf[pI]],
                fCentres[bfI],
                points[bf.prevLabel(pI)]
            );

            vector nPrev = triPrev.normal();
            scalar mPrev = mag(nPrev);

            //- face has zero area
            if( mPrev < VSMALL )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                invertedVertices_.insert(bf[pI]);

                break;
            }
            else
            {
                nPrev /= mPrev;
            }

            //- collocated points
            if( magSqr(triPrev.a() - triPrev.b()) < VSMALL )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                invertedVertices_.insert(bf[pI]);

                break;
            }
            if( magSqr(triPrev.c() - triPrev.a()) < VSMALL )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                invertedVertices_.insert(bf[pI]);

                break;
            }

            //- normal vector is not visible
            if( (nPrev & pNormal) < 0.0 )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                invertedVertices_.insert(bf[pI]);

                break;
            }

            //- check whether the normals of both triangles
            //- point in the same direction
            if( (nNext & nPrev) < 0.0 )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                invertedVertices_.insert(bf[pI]);

                break;
            }
        }
    }

    //- check if there exist concave faces
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        DynList<bool> OkPoints;
        if( !help::isFaceConvexAndOk(bf, points, OkPoints) )
        {
            forAll(OkPoints, pI)
            {
                if( activePointsPtr_ && !(*activePointsPtr_)[bp[bf[pI]]] )
                    continue;

                if( !OkPoints[pI] )
                {
                    # ifdef USE_OMP
                    # pragma omp critical
                    # endif
                    {
                        invertedVertices_.insert(bf[pI]);
                    }
                }
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- exchange global labels of inverted points
        const labelList& bPoints = mse.boundaryPoints();
        const Map<label>& globalToLocal =
            mse.globalToLocalBndPointAddressing();
        const VRWGraph& bpAtProcs = mse.bpAtProcs();
        const DynList<label>& neiProcs = mse.bpNeiProcs();

        std::map<label, labelLongList> shareData;
        forAll(neiProcs, i)
            shareData.insert(std::make_pair(neiProcs[i], labelLongList()));

        forAllConstIter(Map<label>, globalToLocal, iter)
        {
            const label bpI = iter();

            if( !invertedVertices_.found(bPoints[bpI]) )
                continue;

            forAllRow(bpAtProcs, bpI, procI)
            {
                const label neiProc = bpAtProcs(bpI, procI);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                shareData[neiProc].append(iter.key());
            }
        }

        //- exchange data with other processors
        labelLongList receivedData;
        help::exchangeMap(shareData, receivedData);

        forAll(receivedData, i)
        {
            const label bpI = globalToLocal[receivedData[i]];
            invertedVertices_.insert(bPoints[bpI]);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSurfaceCheckInvertedVertices::meshSurfaceCheckInvertedVertices
(
    const meshSurfacePartitioner& mpart
)
:
    surfacePartitioner_(mpart),
    activePointsPtr_(NULL),
    invertedVertices_()
{
    checkVertices();
}

meshSurfaceCheckInvertedVertices::meshSurfaceCheckInvertedVertices
(
    const meshSurfacePartitioner& mpart,
    const boolList& activePoints
)
:
    surfacePartitioner_(mpart),
    activePointsPtr_(&activePoints),
    invertedVertices_()
{
    checkVertices();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceCheckInvertedVertices::~meshSurfaceCheckInvertedVertices()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
