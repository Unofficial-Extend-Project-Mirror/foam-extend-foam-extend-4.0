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

#include "meshSurfaceEngineModifier.H"
#include "polyMeshGenModifier.H"
#include "demandDrivenData.H"

#include "labelledPoint.H"
#include "helperFunctionsPar.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSurfaceEngineModifier::meshSurfaceEngineModifier
(
    meshSurfaceEngine& surfaceEngine
)
:
    surfaceEngine_(surfaceEngine)
{}

meshSurfaceEngineModifier::meshSurfaceEngineModifier
(
    const meshSurfaceEngine& surfaceEngine
)
:
    surfaceEngine_(const_cast<meshSurfaceEngine&>(surfaceEngine))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceEngineModifier::~meshSurfaceEngineModifier()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceEngineModifier::moveBoundaryVertexNoUpdate
(
    const label bpI,
    const point& newP
)
{
    surfaceEngine_.mesh_.points()[surfaceEngine_.boundaryPoints()[bpI]] = newP;
}

void meshSurfaceEngineModifier::moveBoundaryVertex
(
    const label bpI,
    const point& newP
)
{
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    pointFieldPMG& points = surfaceEngine_.mesh_.points();
    points[bPoints[bpI]] = newP;

    if( surfaceEngine_.faceCentresPtr_ )
    {
        vectorField& faceCentres = *surfaceEngine_.faceCentresPtr_;
        const VRWGraph& pFaces = surfaceEngine_.pointFaces();
        const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();

        forAllRow(pFaces, bpI, pfI)
        {
            const label bfI = pFaces(bpI, pfI);

            faceCentres[bfI] = bFaces[bfI].centre(points);
        }
    }

    if( surfaceEngine_.faceNormalsPtr_ )
    {
        vectorField& faceNormals = *surfaceEngine_.faceNormalsPtr_;
        const VRWGraph& pFaces = surfaceEngine_.pointFaces();
        const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();

        forAllRow(pFaces, bpI, pfI)
        {
            const label bfI = pFaces(bpI, pfI);

            faceNormals[bfI] = bFaces[bfI].normal(points);
        }
    }

    if( surfaceEngine_.pointNormalsPtr_ )
    {
        const vectorField& faceNormals = *surfaceEngine_.faceNormalsPtr_;
        const VRWGraph& pFaces = surfaceEngine_.pointFaces();
        const VRWGraph& pPoints = surfaceEngine_.pointPoints();

        vectorField& pn = *surfaceEngine_.pointNormalsPtr_;
        vector n(vector::zero);
        forAllRow(pFaces, bpI, pfI)
            n += faceNormals[pFaces(bpI, pfI)];

        const scalar l = mag(n);
        if( l > VSMALL )
        {
            n /= l;
        }
        else
        {
            n = vector::zero;
        }

        pn[bpI] = n;

        //- change normal of vertices connected to bpI
        forAllRow(pPoints, bpI, ppI)
        {
            const label bpJ = pPoints(bpI, ppI);
            n = vector::zero;
            forAllRow(pFaces, bpJ, pfI)
                n += faceNormals[pFaces(bpJ, pfI)];

            const scalar d = mag(n);
            if( d > VSMALL )
            {
                n /= d;
            }
            else
            {
                n = vector::zero;
            }

            pn[bpJ] = n;
        }
    }
}

void meshSurfaceEngineModifier::syncVerticesAtParallelBoundaries()
{
    if( !Pstream::parRun() )
        return;

    const Map<label>& globalToLocal =
        surfaceEngine_.globalToLocalBndPointAddressing();
    labelLongList syncNodes;
    forAllConstIter(Map<label>, globalToLocal, it)
        syncNodes.append(it());

    syncVerticesAtParallelBoundaries(syncNodes);
}

void meshSurfaceEngineModifier::syncVerticesAtParallelBoundaries
(
    const labelLongList& syncNodes
)
{
    if( !Pstream::parRun() )
        return;

    const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
    const labelList& globalLabel =
        surfaceEngine_.globalBoundaryPointLabel();
    const Map<label>& globalToLocal =
        surfaceEngine_.globalToLocalBndPointAddressing();
    const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.mesh().points();

    std::map<label, LongList<labelledPoint> > exchangeData;
    forAll(neiProcs, i)
        exchangeData.insert
        (
            std::make_pair(neiProcs[i], LongList<labelledPoint>())
        );

    //- construct the map
    forAll(syncNodes, snI)
    {
        const label bpI = syncNodes[snI];

        if( bpAtProcs.sizeOfRow(bpI) == 0 )
            continue;

        point p = points[bPoints[bpI]] / bpAtProcs.sizeOfRow(bpI);
        moveBoundaryVertexNoUpdate(bpI, p);

        forAllRow(bpAtProcs, bpI, i)
        {
            const label neiProc = bpAtProcs(bpI, i);
            if( neiProc == Pstream::myProcNo() )
                continue;

            exchangeData[neiProc].append(labelledPoint(globalLabel[bpI], p));
        }
    }

    //- exchange the data with other processors
    LongList<labelledPoint> receivedData;
    help::exchangeMap(exchangeData, receivedData);

    //- adjust the coordinates
    forAll(receivedData, i)
    {
        const labelledPoint& lp = receivedData[i];
        const label bpI = globalToLocal[lp.pointLabel()];
        const point newP = points[bPoints[bpI]] + lp.coordinates();
        moveBoundaryVertexNoUpdate(bpI, newP);
    }
}

void meshSurfaceEngineModifier::updateGeometry
(
    const labelLongList& updateBndNodes
)
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();
    const VRWGraph& pFaces = surfaceEngine_.pointFaces();
    const labelList& bp = surfaceEngine_.bp();

    boolList updateFaces(bFaces.size(), false);
    # ifdef USE_OMP
    # pragma omp parallel for if( updateBndNodes.size() > 1000 )
    # endif
    forAll(updateBndNodes, i)
    {
        const label bpI = updateBndNodes[i];
        forAllRow(pFaces, bpI, j)
            updateFaces[pFaces(bpI, j)] = true;
    }

    if( surfaceEngine_.faceCentresPtr_ )
    {
        vectorField& faceCentres = *surfaceEngine_.faceCentresPtr_;

        # ifdef USE_OMP
        # pragma omp parallel for if( updateFaces.size() > 1000 ) \
        schedule(dynamic, 100)
        # endif
        forAll(updateFaces, bfI)
        {
            if( updateFaces[bfI] )
                faceCentres[bfI] = bFaces[bfI].centre(points);
        }
    }

    if( surfaceEngine_.faceNormalsPtr_ )
    {
        vectorField& faceNormals = *surfaceEngine_.faceNormalsPtr_;

        # ifdef USE_OMP
        # pragma omp parallel for if( updateFaces.size() > 1000 ) \
        schedule(dynamic, 100)
        # endif
        forAll(updateFaces, bfI)
        {
            if( updateFaces[bfI] )
                faceNormals[bfI] = bFaces[bfI].normal(points);
        }
    }

    if( surfaceEngine_.pointNormalsPtr_ )
    {
        const vectorField& faceNormals = surfaceEngine_.faceNormals();

        boolList updateBndPoint(pFaces.size(), false);
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(updateBndNodes, i)
        {
            const label bpI = updateBndNodes[i];

            forAllRow(pFaces, bpI, pfI)
            {
                const face& bf = bFaces[pFaces(bpI, pfI)];

                forAll(bf, pI)
                    updateBndPoint[bp[bf[pI]]] = true;
            }
        }

        vectorField& pn = *surfaceEngine_.pointNormalsPtr_;
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 100)
        # endif
        forAll(updateBndPoint, bpI)
        {
            if( !updateBndPoint[bpI] )
                continue;

            vector n(vector::zero);
            forAllRow(pFaces, bpI, pfI)
                n += faceNormals[pFaces(bpI, pfI)];

            const scalar l = mag(n);
            if( l > VSMALL )
            {
                n /= l;
            }
            else
            {
                n = vector::zero;
            }

            pn[bpI] = n;
        }

        if( Pstream::parRun() )
        {
            //- update point normals at inter-processor boundaries
            const Map<label>& globalToLocal =
                surfaceEngine_.globalToLocalBndPointAddressing();
            const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
            const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();

            //- make sure that the points ar updated on all processors
            std::map<label, labelLongList> exchangeNodeLabels;
            forAll(neiProcs, i)
                exchangeNodeLabels[neiProcs[i]].clear();

            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                if( updateBndPoint[bpI] )
                {
                    forAllRow(bpAtProcs, bpI, i)
                    {
                        const label neiProc = bpAtProcs(bpI, i);

                        if( neiProc == Pstream::myProcNo() )
                            continue;

                        exchangeNodeLabels[neiProc].append(it.key());
                    }
                }
            }

            labelLongList receivedNodes;
            help::exchangeMap(exchangeNodeLabels, receivedNodes);

            forAll(receivedNodes, i)
                updateBndPoint[globalToLocal[receivedNodes[i]]] = true;


            //- start updating point normals
            std::map<label, LongList<labelledPoint> > exchangeData;
            forAll(neiProcs, i)
                exchangeData[neiProcs[i]].clear();

            //- prepare data for sending
            forAllConstIter(Map<label>, globalToLocal, iter)
            {
                const label bpI = iter();

                if( !updateBndPoint[bpI] )
                    continue;

                vector& n = pn[bpI];
                n = vector::zero;

                forAllRow(pFaces, bpI, pfI)
                    n += faceNormals[pFaces(bpI, pfI)];

                forAllRow(bpAtProcs, bpI, procI)
                {
                    const label neiProc = bpAtProcs(bpI, procI);
                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeData[neiProc].append(labelledPoint(iter.key(), n));
                }
            }

            //- exchange data with other procs
            LongList<labelledPoint> receivedData;
            help::exchangeMap(exchangeData, receivedData);

            forAll(receivedData, i)
            {
                const label bpI = globalToLocal[receivedData[i].pointLabel()];
                pn[bpI] += receivedData[i].coordinates();
            }

            //- normalize vectors
            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                if( !updateBndPoint[bpI] )
                    continue;

                vector normal = pn[bpI];
                const scalar d = mag(normal);
                if( d > VSMALL )
                {
                    normal /= d;
                }
                else
                {
                    normal = vector::zero;
                }

                pn[bpI] = normal;
            }
        }
    }
}

void meshSurfaceEngineModifier::updateGeometry()
{
    labelLongList updateBndNodes(surfaceEngine_.boundaryPoints().size());

    # ifdef USE_OMP
    # pragma omp parallel for if( updateBndNodes.size() > 10000 )
    # endif
    forAll(updateBndNodes, bpI)
        updateBndNodes[bpI] = bpI;

    updateGeometry(updateBndNodes);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
