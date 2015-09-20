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
#include "boundaryLayerOptimisation.H"
#include "meshSurfacePartitioner.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "labelledScalar.H"
#include "refLabelledPoint.H"
#include "refLabelledPointScalar.H"
#include "polyMeshGenAddressing.H"
#include "partTriMesh.H"
#include "partTetMeshSimplex.H"
#include "volumeOptimizer.H"

//#define DEBUGLayer

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayerOptimisation::calculateNormalVectors
(
    const direction eType,
    pointNormalsType& pointPatchNormal
) const
{
    const meshSurfaceEngine& mse = meshSurface();
    const labelList& facePatch = mse.boundaryFacePatches();
    const labelList& bp = mse.bp();
    const VRWGraph& pointFaces = mse.pointFaces();
    const vectorField& fNormals = mse.faceNormals();

    //- calculate point normals with respect to all patches at a point
    pointPatchNormal.clear();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 20)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        if( !(hairEdgeType_[hairEdgeI] & eType) )
            continue;

        const label bpI = bp[hairEdges_[hairEdgeI][0]];

        //- create an entry in a map
        patchNormalType* patchNormalPtr(NULL);
        # ifdef USE_OMP
        # pragma omp critical
            patchNormalPtr = &pointPatchNormal[bpI];
        # else
        patchNormalPtr = &pointPatchNormal[bpI];
        # endif
        patchNormalType& patchNormal = *patchNormalPtr;

        //- sum normals of faces attached to a point
        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            const label patchI = facePatch[bfI];

            if( patchNormal.find(patchI) == patchNormal.end() )
            {
                patchNormal[patchI].first = fNormals[bfI];
                patchNormal[patchI].second = mag(fNormals[bfI]);
            }
            else
            {
                patchNormal[patchI].first += fNormals[bfI];
                patchNormal[patchI].second += mag(fNormals[bfI]);
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- gather information about face normals on other processors
        const Map<label>& globalToLocal =
            mse.globalToLocalBndPointAddressing();
        const DynList<label>& neiProcs = mse.bpNeiProcs();
        const VRWGraph& bpAtProcs = mse.bpAtProcs();

        std::map<label, LongList<refLabelledPointScalar> > exchangeData;
        forAll(neiProcs, i)
            exchangeData[neiProcs[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            if( pointPatchNormal.find(bpI) != pointPatchNormal.end() )
            {
                const patchNormalType& patchNormal = pointPatchNormal[bpI];

                forAllRow(bpAtProcs, bpI, i)
                {
                    const label neiProc = bpAtProcs(bpI, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    forAllConstIter(patchNormalType, patchNormal, pIt)
                        exchangeData[neiProc].append
                        (
                            refLabelledPointScalar
                            (
                                it.key(),
                                labelledPointScalar
                                (
                                    pIt->first,
                                    pIt->second.first,
                                    pIt->second.second
                                )
                            )
                        );
                }
            }
        }

        LongList<refLabelledPointScalar> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const refLabelledPointScalar& rlps = receivedData[i];
            const label bpI = globalToLocal[rlps.objectLabel()];

            patchNormalType& patchNormal = pointPatchNormal[bpI];

            const labelledPointScalar& lps = rlps.lps();
            patchNormal[lps.pointLabel()].first += lps.coordinates();
            patchNormal[lps.pointLabel()].second += lps.scalarValue();
        }
    }

    //- finally, calculate normal vectors
    # ifdef USE_OMP
    # pragma omp parallel
    # pragma omp single nowait
    # endif
    forAllIter(pointNormalsType, pointPatchNormal, it)
    {
        # ifdef USE_OMP
        # pragma omp task firstprivate(it)
        {
        # endif

        patchNormalType& patchNormal = it->second;

        forAllIter(patchNormalType, patchNormal, pIt)
        {
            pIt->second.first /= pIt->second.second;
            //pIt->second.first /= (mag(pIt->second.first) + VSMALL);
        }

        # ifdef USE_OMP
        }
        # endif
    }
}

void boundaryLayerOptimisation::calculateNormalVectorsSmother
(
    const direction eType,
    pointNormalsType& pointPatchNormal
)
{
    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const meshSurfaceEngine& mse = mPart.surfaceEngine();
    const pointFieldPMG& points = mse.points();
    const labelList& bp = mse.bp();

    partTriMesh triMesh(mPart);

    const pointField& triMeshPoints = triMesh.points();
    const VRWGraph& pTriangles = triMesh.pointTriangles();
    const LongList<labelledTri>& triangles = triMesh.triangles();
    const labelList& triPointLabel = triMesh.meshSurfacePointLabelInTriMesh();
    const labelLongList& surfPointLabel = triMesh.pointLabelInMeshSurface();

    Info << "Calculating normals using smoother " << endl;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(hairEdgeType_, heI)
    {
        if( !(hairEdgeType_[heI] & eType) )
            continue;

        const edge& he = hairEdges_[heI];
        const label bpI = bp[he.start()];

        const label triPointI = triPointLabel[bpI];

        //- create an entry in a map
        patchNormalType* patchNormalPtr(NULL);
        # ifdef USE_OMP
        # pragma omp critical
            patchNormalPtr = &pointPatchNormal[bpI];
        # else
        patchNormalPtr = &pointPatchNormal[bpI];
        # endif

        patchNormalType& patchNormal = *patchNormalPtr;

        //- find patches at this point
        DynList<label> patchesAtPoint;
        forAllRow(pTriangles, triPointI, ptI)
        {
            patchesAtPoint.appendIfNotIn
            (
                triangles[pTriangles(triPointI, ptI)].region()
            );
        }

        forAll(patchesAtPoint, ptchI)
        {
            const label patchI = patchesAtPoint[ptchI];

            DynList<point, 128> pts(2);
            DynList<partTet, 128> tets;

            //- create points
            pts[0] = points[he.start()];
            pts[1] = points[he.end()];

            Map<label> bpToSimplex;
            bpToSimplex.insert(bpI, 0);

            forAllRow(pTriangles, triPointI, ptI)
            {
                const labelledTri& tri = triangles[pTriangles(triPointI, ptI)];

                if( tri.region() == patchI )
                {
                    //- create points originating from triangles
                    FixedList<label, 3> triLabels;
                    forAll(tri, pI)
                    {
                        const label spLabel = tri[pI];
                        const label bpLabel = surfPointLabel[spLabel];

                        if( bpLabel < 0 )
                        {
                            triLabels[pI] = pts.size();
                            pts.append(triMeshPoints[spLabel]);
                            continue;
                        }

                        if( !bpToSimplex.found(bpLabel) )
                        {
                            bpToSimplex.insert(bpLabel, pts.size());
                            pts.append(triMeshPoints[spLabel]);
                        }

                        triLabels[pI] = bpToSimplex[bpLabel];
                    }

                    //- create a new tetrahedron
                    tets.append
                    (
                        partTet(triLabels[2], triLabels[1], triLabels[0], 1)
                    );
                }
            }

            //- create partTeMeshSimplex
            partTetMeshSimplex simplex(pts, tets, 1);

            //- activate volume optimizer
            volumeOptimizer vOpt(simplex);

            vOpt.optimizeNodePosition();

            const point newP = simplex.centrePoint();

            vector n = -1.0 * (newP - pts[0]);
            const scalar magN = (mag(n) + VSMALL);

            patchNormal[patchI].first = (n / magN);
            patchNormal[patchI].second = magN;
        }
    }

    Info << "Finished calculating normals using smoother " << endl;
}

void boundaryLayerOptimisation::calculateHairVectorsAtTheBoundary
(
    vectorField& hairVecs
)
{
    //- set the size of hairVecs
    hairVecs.setSize(hairEdges_.size());

    const pointFieldPMG& points = mesh_.points();

    const meshSurfaceEngine& mse = meshSurface();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelList& bp = mse.bp();
    const VRWGraph& bpEdges = mse.boundaryPointEdges();
    const edgeList& edges = mse.edges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        const direction hairType = hairEdgeType_[hairEdgeI];

        if( hairType & BOUNDARY )
        {
            const edge& he = hairEdges_[hairEdgeI];
            vector& hv = hairVecs[hairEdgeI];
            hv = vector::zero;

            if( hairType & FEATUREEDGE )
            {
                //- do not modify hair vectors at feature edges
                hv = he.vec(points);
            }
            else if( hairType & (ATEDGE|ATCORNER) )
            {
                //- this is a case of O-layer at a corner or feature edge

                //- find the surface edges corresponding to the hair edge
                label beI(-1);

                const label bps = bp[he.start()];
                forAllRow(bpEdges, bps, bpeI)
                {
                    const label beJ = bpEdges(bps, bpeI);

                    if( edges[beJ] == he )
                    {
                        beI = beJ;
                        continue;
                    }
                }

                if( beI < 0 )
                    FatalErrorIn
                    (
                        "boundaryLayerOptimisation::"
                        "calculateHairVectorsAtTheBoundary(vectorField&)"
                    ) << "Cannot find hair edge "
                      << hairEdgeI << abort(FatalError);

                //- find the vector at the same angle from both feature edges
                forAllRow(edgeFaces, beI, befI)
                {
                    const face& bf = bFaces[edgeFaces(beI, befI)];
                    const vector fNormal = bf.normal(points);

                    const label pos = bf.which(he.start());

                    if( pos < 0 )
                        FatalErrorIn
                        (
                            "boundaryLayerOptimisation::"
                            "calculateHairVectorsAtTheBoundary(vectorField&)"
                        ) << "Cannot find hair edge "
                          << hairEdgeI << " in face " << bf
                          << abort(FatalError);

                    if( he.end() == bf.prevLabel(pos) )
                    {
                        const edge fe = bf.faceEdge(pos);
                        const vector ev = fe.vec(points);

                        vector hev = fNormal ^ ev;
                        hev /= (mag(hev) + VSMALL);

                        hv += hev;
                    }
                    else
                    {
                        const edge fe = bf.faceEdge(bf.rcIndex(pos));
                        const vector ev = fe.vec(points);

                        vector hev = fNormal ^ ev;
                        hev /= (mag(hev) + VSMALL);

                        hv += hev;
                    }
                }
            }
            else
            {
                FatalErrorIn
                (
                    "boundaryLayerOptimisation::"
                    "calculateHairVectorsAtTheBoundary(vectorField&)"
                ) << "Invalid hair type " << label(hairType)
                  << abort(FatalError);
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- collect data at inter-processor boundaries
        const Map<label>& globalToLocal =
            mse.globalToLocalBndPointAddressing();

        const edgeList& edges = mesh_.addressingData().edges();
        const VRWGraph& pointEdges = mesh_.addressingData().pointEdges();
        const VRWGraph& edgesAtProcs =
            mesh_.addressingData().edgeAtProcs();
        const labelLongList& globalEdgeLabel =
            mesh_.addressingData().globalEdgeLabel();
        const Map<label>& globalToLocalEdge =
            mesh_.addressingData().globalToLocalEdgeAddressing();
        const DynList<label>& eNeiProcs =
            mesh_.addressingData().edgeNeiProcs();

        std::map<label, LongList<labelledPoint> > exchangeData;
        forAll(eNeiProcs, i)
            exchangeData[eNeiProcs[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            forAllRow(hairEdgesAtBndPoint_, bpI, i)
            {
                const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);
                const edge& he = hairEdges_[hairEdgeI];

                const direction eType = hairEdgeType_[hairEdgeI];

                //- filte out edges which are not relevant
                if( !(eType & BOUNDARY) )
                    continue;

                forAllRow(pointEdges, he.start(), peI)
                {
                    const label edgeI = pointEdges(he.start(), peI);

                    if( he == edges[edgeI] )
                    {
                        labelledPoint lp
                        (
                            globalEdgeLabel[edgeI],
                            hairVecs[hairEdgeI]
                        );

                        forAllRow(edgesAtProcs, edgeI, j)
                        {
                            const label neiProc = edgesAtProcs(edgeI, j);

                            if( neiProc == Pstream::myProcNo() )
                                continue;

                            LongList<labelledPoint>& dts =
                                exchangeData[neiProc];

                            dts.append(lp);
                        }
                    }
                }
            }
        }

        LongList<labelledPoint> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const labelledPoint& lp = receivedData[i];
            const label edgeI = globalToLocalEdge[lp.pointLabel()];
            const edge& e = edges[edgeI];

            bool found(false);
            for(label pI=0;pI<2;++pI)
            {
                const label bpI = bp[e[pI]];

                if( bpI < 0 )
                    continue;

                forAllRow(hairEdgesAtBndPoint_, bpI, i)
                {
                    const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);

                    if( hairEdges_[hairEdgeI] == e )
                    {
                        hairVecs[hairEdgeI] += lp.coordinates();

                        found = true;
                        break;
                    }
                }

                if( found )
                    break;
            }
        }
    }

    //- calculate new normal vectors
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(hairVecs, hairEdgeI)
    {
        if( hairEdgeType_[hairEdgeI] & BOUNDARY )
            hairVecs[hairEdgeI] /= (mag(hairVecs[hairEdgeI]) + VSMALL);
    }

    # ifdef DEBUGLayer
    Info << "Saving bnd hair vectors" << endl;
    writeHairEdges("bndHairVectors.vtk", (BOUNDARY | ATEDGE), hairVecs);
    # endif
}

void boundaryLayerOptimisation::optimiseHairNormalsAtTheBoundary()
{
    pointFieldPMG& points = mesh_.points();

    //- calculate direction of hair vector based on the surface normal
    const meshSurfaceEngine& mse = meshSurface();
    const labelList& bp = mse.bp();
    const VRWGraph& pFaces = mse.pointFaces();
    const faceList::subList& bFaces = mse.boundaryFaces();

    //- calculate hair vectors
    //- they point in the normal direction to the surface
    vectorField hairVecs(hairEdges_.size());

    if( reCalculateNormals_ )
    {
        //- calulate new normal vectors
        calculateHairVectorsAtTheBoundary(hairVecs);
    }
    else
    {
        if( nSmoothNormals_ == 0 )
            return;

        //- keep existing hair vectors
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdges_, heI)
            hairVecs[heI] = hairEdges_[heI].vec(points);
    }

    Info << "Smoothing boundary hair vectors" << endl;

    //- smooth the variation of normals to reduce the twisting of faces
    label nIter(0);

    while( nIter++ < nSmoothNormals_ )
    {
        vectorField newNormals(hairVecs.size());

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdgesNearHairEdge_, hairEdgeI)
        {
            vector& newNormal = newNormals[hairEdgeI];
            newNormal = vector::zero;

            const direction eType = hairEdgeType_[hairEdgeI];

            if( eType & BOUNDARY )
            {
                if( eType & (FEATUREEDGE | ATCORNER) )
                {
                    //- hair vectors at feature edges must not be modified
                    newNormal += hairVecs[hairEdgeI];
                }
                else if( eType & ATEDGE )
                {
                    const edge& he = hairEdges_[hairEdgeI];

                    DynList<label, 2> edgeFaces;
                    const label bps = bp[he.start()];
                    forAllRow(pFaces, bps, pfI)
                    {
                        const label bfI = pFaces(bps, pfI);
                        const face& bf = bFaces[bfI];

                        forAll(bf, eI)
                        {
                            if( bf.faceEdge(eI) == he )
                            {
                                edgeFaces.append(bfI);
                            }
                        }
                    }

                    //- find the best fitting vector
                    //- at the surface of the mesh
                    forAllRow(hairEdgesNearHairEdge_, hairEdgeI, nheI)
                    {
                        const label hairEdgeJ =
                            hairEdgesNearHairEdge_(hairEdgeI, nheI);

                        //- check if the neighbour hair edge shares a boundary
                        //- face with the current hair edge
                        bool useNeighbour(false);
                        const edge& nhe = hairEdges_[hairEdgeJ];
                        forAll(edgeFaces, efI)
                        {
                            const face& bf = bFaces[edgeFaces[efI]];

                            forAll(bf, eI)
                            {
                                if( bf.faceEdge(eI) == nhe )
                                {
                                    useNeighbour = true;
                                    break;
                                }
                            }
                        }

                        if( useNeighbour )
                            newNormal += hairVecs[hairEdgeJ];
                    }
                }
                else
                {
                    FatalErrorIn
                    (
                        "void boundaryLayerOptimisation::"
                        "optimiseHairNormalsAtTheBoundary(const label)"
                    ) << "Cannot smooth hair with type " << label(eType)
                      << abort(FatalError);
                }
            }
        }

        if( Pstream::parRun() )
        {
            //- collect data at inter-processor boundaries
            const Map<label>& globalToLocal =
                mse.globalToLocalBndPointAddressing();

            const edgeList& edges = mesh_.addressingData().edges();
            const VRWGraph& pointEdges = mesh_.addressingData().pointEdges();
            const VRWGraph& edgesAtProcs =
                mesh_.addressingData().edgeAtProcs();
            const labelLongList& globalEdgeLabel =
                mesh_.addressingData().globalEdgeLabel();
            const Map<label>& globalToLocalEdge =
                mesh_.addressingData().globalToLocalEdgeAddressing();
            const DynList<label>& eNeiProcs =
                mesh_.addressingData().edgeNeiProcs();

            std::map<label, LongList<labelledPoint> > exchangeData;
            forAll(eNeiProcs, i)
                exchangeData[eNeiProcs[i]].clear();

            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                forAllRow(hairEdgesAtBndPoint_, bpI, i)
                {
                    const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);
                    const edge& he = hairEdges_[hairEdgeI];

                    const direction eType = hairEdgeType_[hairEdgeI];

                    //- filte out edges which are not relevant
                    if( !(eType & BOUNDARY) )
                        continue;
                    if( !(eType & ATEDGE) )
                        continue;
                    if( eType & (FEATUREEDGE|ATCORNER) )
                        continue;

                    forAllRow(pointEdges, he.start(), peI)
                    {
                        const label edgeI = pointEdges(he.start(), peI);

                        if( he == edges[edgeI] )
                        {
                            labelledPoint lp
                            (
                                globalEdgeLabel[edgeI],
                                hairVecs[hairEdgeI]
                            );

                            forAllRow(edgesAtProcs, edgeI, j)
                            {
                                const label neiProc = edgesAtProcs(edgeI, j);

                                if( neiProc == Pstream::myProcNo() )
                                    continue;

                                LongList<labelledPoint>& dts =
                                    exchangeData[neiProc];

                                dts.append(lp);
                            }
                        }
                    }
                }
            }

            LongList<labelledPoint> receivedData;
            help::exchangeMap(exchangeData, receivedData);

            forAll(receivedData, i)
            {
                const labelledPoint& lp = receivedData[i];
                const label edgeI = globalToLocalEdge[lp.pointLabel()];
                const edge& e = edges[edgeI];

                bool found(false);
                for(label pI=0;pI<2;++pI)
                {
                    const label bpI = bp[e[pI]];

                    if( bpI < 0 )
                        continue;

                    forAllRow(hairEdgesAtBndPoint_, bpI, i)
                    {
                        const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);

                        if( hairEdges_[hairEdgeI] == e )
                        {
                            hairVecs[hairEdgeI] += lp.coordinates();

                            found = true;
                            break;
                        }
                    }

                    if( found )
                        break;
                }
            }
        }

        //- calculate new normal vectors
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 100)
        # endif
        forAll(newNormals, heI)
        {
            if( hairEdgeType_[heI] & BOUNDARY )
            {
                newNormals[heI] /= (mag(newNormals[heI]) + VSMALL);
                newNormals[heI] = 0.5 * (newNormals[heI] + hairVecs[heI]);
                newNormals[heI] /= (mag(newNormals[heI]) + VSMALL);
            }
        }

        //- transfer new hair vectors to the hairVecs list
        hairVecs.transfer(newNormals);

        # ifdef DEBUGLayer
        if( true )
        {
            writeHairEdges
            (
                "bndHairVectors_"+help::scalarToText(nIter)+".vtk",
                (BOUNDARY | ATEDGE),
                hairVecs
            );
        }
        # endif
    }

    //- move vertices to the new locations
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        if( hairEdgeType_[hairEdgeI] & BOUNDARY )
        {
            const edge& he = hairEdges_[hairEdgeI];

            const vector& hv = hairVecs[hairEdgeI];

            points[he.end()] = points[he.start()] + hv * he.mag(points);
        }
    }

    Info << "Finished smoothing boundary hair vectors" << endl;
}

void boundaryLayerOptimisation::optimiseHairNormalsInside()
{
    pointFieldPMG& points = mesh_.points();

    //- calculate direction of hair vector based on the surface normal
    const meshSurfaceEngine& mse = meshSurface();
    const labelList& bp = mse.bp();
    const VRWGraph& bpEdges = mse.boundaryPointEdges();
    const edgeList& edges = mse.edges();

    //- they point in the normal direction to the surface
    vectorField hairVecs(hairEdges_.size());

    if( reCalculateNormals_ )
    {
        //- calculate point normals with respect to all patches at a point
        pointNormalsType pointPatchNormal;
        calculateNormalVectors(INSIDE, pointPatchNormal);
        //calculateNormalVectorsSmother(INSIDE, pointPatchNormal);

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdges_, hairEdgeI)
        {
            const direction hairType = hairEdgeType_[hairEdgeI];

            if( hairType & INSIDE )
            {
                vector& hv = hairVecs[hairEdgeI];
                hv = vector::zero;

                const label bpI = bp[hairEdges_[hairEdgeI].start()];

                label counter(0);
                const patchNormalType& patchNormals = pointPatchNormal[bpI];
                forAllConstIter(patchNormalType, patchNormals, pIt)
                {
                    hv -= pIt->second.first;
                    ++counter;
                }

                if( counter == 0 )
                {
                    FatalErrorIn
                    (
                        "void boundaryLayerOptimisation::"
                        "optimiseHairNormalsInside()"
                    ) << "No valid patches for boundary point "
                      << bp[hairEdges_[hairEdgeI].start()] << abort(FatalError);
                }

                hv /= (mag(hv) + VSMALL);
            }
            else
            {
                //- initialise boundary hair vectors. They influence internal
                //- hairs connected to them
                vector hvec = hairEdges_[hairEdgeI].vec(points);
                hvec /= (mag(hvec) + VSMALL);
                hairVecs[hairEdgeI] = hvec;
            }
        }
    }
    else
    {
        if( nSmoothNormals_ == 0 )
            return;

        Info << "Using existing hair vectors" << endl;

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdges_, hairEdgeI)
        {
            //- initialise boundary hair vectors.
            vector hvec = hairEdges_[hairEdgeI].vec(points);
            hvec /= (mag(hvec) + VSMALL);
            hairVecs[hairEdgeI] = hvec;
        }
    }

    # ifdef DEBUGLayer
    writeHairEdges("insideHairVectors.vtk", (INSIDE|BOUNDARY), hairVecs);
    # endif

    Info << "Smoothing internal hair vectors" << endl;

    //- smooth the variation of normals to reduce twisting of faces
    label nIter(0);

    while( nIter++ < nSmoothNormals_ )
    {
        vectorField newNormals(hairVecs.size());

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdgesNearHairEdge_, hairEdgeI)
        {
            vector& newNormal = newNormals[hairEdgeI];
            newNormal = vector::zero;

            const direction eType = hairEdgeType_[hairEdgeI];

            const edge& he = hairEdges_[hairEdgeI];
            const vector& heVec = hairVecs[hairEdgeI];

            if( eType & INSIDE )
            {
                if( eType & ATCORNER )
                {
                    //- hair vectors at feature edges must not be modified
                    newNormal = hairVecs[hairEdgeI];
                }
                else if( eType & ATEDGE )
                {
                    forAllRow(hairEdgesNearHairEdge_, hairEdgeI, nheI)
                    {
                        const label hairEdgeJ =
                            hairEdgesNearHairEdge_(hairEdgeI, nheI);

                        const edge& nhe = hairEdges_[hairEdgeJ];
                        const vector& nhVec = hairVecs[hairEdgeJ];

                        vector n = nhVec ^ (points[nhe[0]] - points[he[0]]);
                        n /= (mag(n) + VSMALL);

                        vector newVec = heVec - (heVec & n) * n;
                        newVec /= (mag(newVec) + VSMALL);

                        scalar weight = 1.0;

                        if( Pstream::parRun() )
                        {
                            //- edges at inter-processor boundaries contribute
                            //- at two sides are given weight 0.5
                            const edge be(he[0], nhe[0]);
                            const label bpI = bp[he[0]];

                            forAllRow(bpEdges, bpI, bpeI)
                            {
                                const edge& bndEdge = edges[bpEdges(bpI, bpeI)];

                                if( bndEdge == be )
                                {
                                    weight = 0.5;
                                    break;
                                }
                            }
                        }

                        newNormal += weight * newVec;
                    }
                }
                else
                {
                    //- find the best fitting vector
                    //- at the surface of the mesh
                    forAllRow(hairEdgesNearHairEdge_, hairEdgeI, nheI)
                    {
                        const label hairEdgeJ =
                            hairEdgesNearHairEdge_(hairEdgeI, nheI);

                        scalar weight = 1.0;

                        if( Pstream::parRun() )
                        {
                            //- edges at inter-processor boundaries contribute
                            //- at two sides are given weight 0.5
                            const edge& nhe = hairEdges_[hairEdgeJ];
                            const edge be(he[0], nhe[0]);
                            const label bpI = bp[he[0]];

                            forAllRow(bpEdges, bpI, bpeI)
                            {
                                const edge& bndEdge = edges[bpEdges(bpI, bpeI)];

                                if( bndEdge == be )
                                {
                                    weight = 0.5;
                                    break;
                                }
                            }
                        }

                        newNormal += weight * hairVecs[hairEdgeJ];
                    }
                }
            }
            else
            {
                //- copy the existing hair vector
                newNormal = hairVecs[hairEdgeI];
            }
        }

        if( Pstream::parRun() )
        {
            //- collect data at inter-processor boundaries
            const Map<label>& globalToLocal =
                mse.globalToLocalBndPointAddressing();

            const edgeList& edges = mesh_.addressingData().edges();
            const VRWGraph& pointEdges = mesh_.addressingData().pointEdges();
            const VRWGraph& edgesAtProcs =
                mesh_.addressingData().edgeAtProcs();
            const labelLongList& globalEdgeLabel =
                mesh_.addressingData().globalEdgeLabel();
            const Map<label>& globalToLocalEdge =
                mesh_.addressingData().globalToLocalEdgeAddressing();
            const DynList<label>& eNeiProcs =
                mesh_.addressingData().edgeNeiProcs();

            std::map<label, LongList<labelledPoint> > exchangeData;
            forAll(eNeiProcs, i)
                exchangeData[eNeiProcs[i]].clear();

            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                forAllRow(hairEdgesAtBndPoint_, bpI, i)
                {
                    const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);
                    const edge& he = hairEdges_[hairEdgeI];

                    const direction eType = hairEdgeType_[hairEdgeI];

                    //- handle boundary points, only
                    if( !(eType & INSIDE) || (eType & ATCORNER) )
                        continue;

                    forAllRow(pointEdges, he.start(), peI)
                    {
                        const label edgeI = pointEdges(he.start(), peI);

                        if( he == edges[edgeI] )
                        {
                            labelledPoint lp
                            (
                                globalEdgeLabel[edgeI],
                                hairVecs[hairEdgeI]
                            );

                            forAllRow(edgesAtProcs, edgeI, j)
                            {
                                const label neiProc = edgesAtProcs(edgeI, j);

                                if( neiProc == Pstream::myProcNo() )
                                    continue;

                                LongList<labelledPoint>& dts =
                                    exchangeData[neiProc];

                                dts.append(lp);
                            }
                        }
                    }
                }
            }

            LongList<labelledPoint> receivedData;
            help::exchangeMap(exchangeData, receivedData);

            forAll(receivedData, i)
            {
                const labelledPoint& lp = receivedData[i];
                const label edgeI = globalToLocalEdge[lp.pointLabel()];
                const edge& e = edges[edgeI];

                bool found(false);
                for(label pI=0;pI<2;++pI)
                {
                    const label bpI = bp[e[pI]];

                    if( bpI < 0 )
                        continue;

                    forAllRow(hairEdgesAtBndPoint_, bpI, i)
                    {
                        const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);

                        if( hairEdges_[hairEdgeI] == e )
                        {
                            hairVecs[hairEdgeI] += lp.coordinates();

                            found = true;
                            break;
                        }
                    }

                    if( found )
                        break;
                }
            }
        }

        //- calculate new normal vectors
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 100)
        # endif
        forAll(newNormals, hairEdgeI)
        {
            if( hairEdgeType_[hairEdgeI] & INSIDE )
                newNormals[hairEdgeI] /= (mag(newNormals[hairEdgeI]) + VSMALL);
        }

        //- transfer new hair vectors to the hairVecs list
        hairVecs.transfer(newNormals);

        # ifdef DEBUGLayer
        if( true )
        {
            writeHairEdges
            (
                "insideHairVectors_"+help::scalarToText(nIter)+".vtk",
                INSIDE,
                hairVecs
            );
        }
        # endif
    }

    //- move vertices to the new locations
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        if( hairEdgeType_[hairEdgeI] & INSIDE )
        {
            const edge& he = hairEdges_[hairEdgeI];

            const vector& hv = hairVecs[hairEdgeI];

            points[he.end()] = points[he.start()] + hv * he.mag(points);
        }
    }

    Info << "Finished smoothing internal hair vectors" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
