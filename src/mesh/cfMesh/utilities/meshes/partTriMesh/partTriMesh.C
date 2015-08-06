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
#include "partTriMesh.H"
#include "meshSurfacePartitioner.H"
#include "VRWGraphList.H"
#include "triSurfModifier.H"
#include "helperFunctions.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

partTriMesh::partTriMesh(const meshSurfacePartitioner& mPart)
:
    mPart_(mPart),
    surf_(),
    pointLabelInMeshSurface_(),
    meshSurfacePointLabelInTriMesh_(),
    pointType_(),
    globalPointLabelPtr_(NULL),
    pAtProcsPtr_(NULL),
    globalToLocalPointAddressingPtr_(NULL),
    neiProcsPtr_(NULL),
    pAtParallelBoundariesPtr_(NULL),
    pAtBufferLayersPtr_(NULL)
{
    const meshSurfaceEngine& meshSurface = mPart.surfaceEngine();
    List<direction> useFace(meshSurface.boundaryFaces().size(), direction(1));

    createPointsAndTrias(useFace);
}

partTriMesh::partTriMesh
(
    const meshSurfacePartitioner& mPart,
    const labelHashSet& invertedPoints,
    const label additionalLayers
)
:
    mPart_(mPart),
    surf_(),
    pointLabelInMeshSurface_(),
    meshSurfacePointLabelInTriMesh_(),
    pointType_(),
    globalPointLabelPtr_(NULL),
    pAtProcsPtr_(NULL),
    globalToLocalPointAddressingPtr_(NULL),
    neiProcsPtr_(NULL),
    pAtParallelBoundariesPtr_(NULL),
    pAtBufferLayersPtr_(NULL)
{
    const meshSurfaceEngine& meshSurface = mPart.surfaceEngine();
    const VRWGraph& pointFaces = meshSurface.pointFaces();
    const faceList::subList& bFaces = meshSurface.boundaryFaces();
    const labelList& bPoints = meshSurface.boundaryPoints();
    const labelList& bp = meshSurface.bp();

    List<direction> useFace(bFaces.size(), direction(0));

    //- select cells containing at least one vertex of the bad faces
    forAll(pointFaces, bpI)
        if( invertedPoints.found(bPoints[bpI]) )
        {
            forAllRow(pointFaces, bpI, pfI)
                useFace[pointFaces(bpI, pfI)] = 1;
        }

    //- add additional layer of cells
    for(direction layerI=1;layerI<(additionalLayers+direction(1));++layerI)
    {
        forAll(useFace, bfI)
            if( useFace[bfI] == layerI )
            {
                const face& bf = bFaces[bfI];

                forAll(bf, pI)
                {
                    const label bpI = bp[bf[pI]];

                    forAllRow(pointFaces, bpI, pfI)
                    {
                        const label fLabel = pointFaces(bpI, pfI);

                        if( !useFace[fLabel] )
                            useFace[fLabel] = layerI + 1;
                    }
                }
            }

        if( Pstream::parRun() )
        {
            const labelList& globalPointLabel =
                meshSurface.globalBoundaryPointLabel();
            const VRWGraph& pProcs = meshSurface.bpAtProcs();
            const Map<label>& globalToLocal =
                meshSurface.globalToLocalBndPointAddressing();

            std::map<label, labelLongList> eData;
            forAllConstIter(Map<label>, globalToLocal, iter)
            {
                const label bpI = iter();

                forAllRow(pProcs, bpI, procI)
                {
                    const label neiProc = pProcs(bpI, procI);
                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    if( eData.find(neiProc) == eData.end() )
                    {
                        eData.insert
                        (
                            std::make_pair(neiProc, labelLongList())
                        );
                    }

                    forAllRow(pointFaces, bpI, pfI)
                        if( useFace[pointFaces(bpI, pfI)] == layerI )
                        {
                            eData[neiProc].append(globalPointLabel[bpI]);
                            break;
                        }
                }
            }

            //- exchange data with other processors
            labelLongList receivedData;
            help::exchangeMap(eData, receivedData);

            forAll(receivedData, i)
            {
                const label bpI = globalToLocal[receivedData[i]];

                forAllRow(pointFaces, bpI, pfI)
                {
                    const label fLabel = pointFaces(bpI, pfI);
                    if( !useFace[fLabel] )
                        useFace[fLabel] = layerI + 1;
                }
            }
        }
    }

    createPointsAndTrias(useFace);
}

partTriMesh::~partTriMesh()
{
    deleteDemandDrivenData(globalPointLabelPtr_);
    deleteDemandDrivenData(pAtProcsPtr_);
    deleteDemandDrivenData(globalToLocalPointAddressingPtr_);
    deleteDemandDrivenData(neiProcsPtr_);
    deleteDemandDrivenData(pAtParallelBoundariesPtr_);
    deleteDemandDrivenData(pAtBufferLayersPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void partTriMesh::updateVertex(const label pI, const point& newP)
{
    triSurfModifier sMod(surf_);
    pointField& pts = sMod.pointsAccess();
    const VRWGraph& pointFacets = surf_.pointFacets();

    pts[pI] = newP;

    if( pointType_[pI] & FACECENTRE )
    {
        Warning << "Smoothing auxiliary vertex."
            << " This has no effect on the original mesh" << endl;
        return;
    }

    //- find face centres attached
    DynList<label> helper;
    forAllRow(pointFacets, pI, ptI)
    {
        const label centreI = surf_[pointFacets(pI, ptI)][2];
        if( pointType_[centreI] & FACECENTRE )
            helper.appendIfNotIn(centreI);
    }

    //- update coordinates of FACECENTRE vertices
    forAll(helper, i)
    {
        const label centreI = helper[i];

        point centre(vector::zero);
        scalar faceArea(0.0);
        forAllRow(pointFacets, centreI, ptI)
        {
            const labelledTri& tri = surf_[pointFacets(centreI, ptI)];
            point c(vector::zero);
            for(label i=0;i<3;++i)
                c += pts[tri[i]];
            c /= 3;
            const scalar area = tri.mag(pts) + VSMALL;

            centre += c * area;
            faceArea += area;
        }

        pts[centreI] = centre / faceArea;
    }
}

void partTriMesh::updateVerticesSMP(const List<LongList<labelledPoint> >& np)
{
    triSurfModifier sMod(surf_);
    pointField& pts = sMod.pointsAccess();
    const VRWGraph& pointFacets = surf_.pointFacets();

    List<direction> updateType(pts.size(), direction(0));

    # ifdef USE_OMP
    # pragma omp parallel num_threads(np.size())
    # endif
    {
        # ifdef USE_OMP
        const LongList<labelledPoint>& newPoints = np[omp_get_thread_num()];
        # else
        const LongList<labelledPoint>& newPoints = np[0];
        # endif

        forAll(newPoints, i)
        {
            const labelledPoint& lp = newPoints[i];
            const label pointI = lp.pointLabel();

            pts[pointI] = lp.coordinates();
            updateType[pointI] |= SMOOTH;

            forAllRow(pointFacets, pointI, ptI)
            {
                const labelledTri& tri = surf_[pointFacets(pointI, ptI)];

                if( pointType_[tri[2]] & FACECENTRE )
                    updateType[tri[2]] |= FACECENTRE;
            }
        }
    }

    //- update coordinates of buffer layer points
    if( Pstream::parRun() )
    {
        const labelLongList& bufferLayerPoints = this->bufferLayerPoints();
        const VRWGraph& pProcs = this->pointAtProcs();
        const labelLongList& globalPointLabel = this->globalPointLabel();
        const Map<label>& globalToLocal = this->globalToLocalPointAddressing();
        const DynList<label>& neiProcs = this->neiProcs();

        //- create the map
        std::map<label, LongList<labelledPoint> > exchangeData;
        forAll(neiProcs, i)
            exchangeData.insert
            (
                std::make_pair(neiProcs[i], LongList<labelledPoint>())
            );

        //- add points into the map
        forAll(bufferLayerPoints, pI)
        {
            const label pointI = bufferLayerPoints[pI];

            if( !(updateType[pointI] & SMOOTH) )
                continue;

            forAllRow(pProcs, pointI, i)
            {
                const label neiProc = pProcs(pointI, i);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                exchangeData[neiProc].append
                (
                    labelledPoint(globalPointLabel[pointI], pts[pointI])
                );
            }
        }

        LongList<labelledPoint> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const labelledPoint& lp = receivedData[i];

            const label triPointI = globalToLocal[lp.pointLabel()];
            pts[triPointI] = lp.coordinates();

            forAllRow(pointFacets, triPointI, ptI)
            {
                const labelledTri& tri = surf_[pointFacets(triPointI, ptI)];

                if( pointType_[tri[2]] & FACECENTRE )
                    updateType[tri[2]] |= FACECENTRE;
            }
        }
    }

    //- update coordinates of FACECENTRE vertices
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 20)
    # endif
    forAll(updateType, pI)
    {
        if( updateType[pI] & FACECENTRE )
        {
            point centre(vector::zero);
            scalar faceArea(0.0);
            forAllRow(pointFacets, pI, ptI)
            {
                const labelledTri& tri = surf_[pointFacets(pI, ptI)];
                point c(vector::zero);
                for(label i=0;i<3;++i)
                    c += pts[tri[i]];
                c /= 3;
                const scalar area = tri.mag(pts) + VSMALL;

                centre += c * area;
                faceArea += area;
            }

            pts[pI] = centre / faceArea;
        }
    }

    if( Pstream::parRun() )
        updateBufferLayers();
}

void partTriMesh::updateVertices()
{
    const meshSurfaceEngine& mse = mPart_.surfaceEngine();
    const labelList& bPoints = mse.boundaryPoints();
    labelLongList movedPoints(bPoints.size());
    forAll(movedPoints, i)
        movedPoints[i] = i;

    updateVertices(movedPoints);
}

void partTriMesh::updateVertices(const labelLongList& movedPoints)
{
    const meshSurfaceEngine& mse = mPart_.surfaceEngine();
    const pointFieldPMG& points = mse.points();
    const labelList& bPoints = mse.boundaryPoints();

    triSurfModifier sMod(surf_);
    pointField& pts = sMod.pointsAccess();
    const VRWGraph& pointFacets = surf_.pointFacets();

    List<direction> updateType(pts.size(), direction(0));

    //- update coordinates of vertices which exist in the surface
    //- of the volume mesh
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(movedPoints, i)
    {
        const label bpI = movedPoints[i];
        const label pointI = bPoints[bpI];
        const label triPointI = meshSurfacePointLabelInTriMesh_[bpI];

        if( triPointI < 0 )
            continue;

        pts[triPointI] = points[pointI];
        updateType[triPointI] |= SMOOTH;

        forAllRow(pointFacets, triPointI, ptI)
        {
            const labelledTri& tri = surf_[pointFacets(triPointI, ptI)];

            if( pointType_[tri[2]] & FACECENTRE )
                updateType[tri[2]] |= FACECENTRE;
        }
    }

    //- update coordinates of buffer layer points
    if( Pstream::parRun() )
    {
        const labelLongList& bufferLayerPoints = this->bufferLayerPoints();
        const VRWGraph& pProcs = this->pointAtProcs();
        const labelLongList& globalPointLabel = this->globalPointLabel();
        const Map<label>& globalToLocal = this->globalToLocalPointAddressing();
        const DynList<label>& neiProcs = this->neiProcs();

        //- create the map
        std::map<label, LongList<labelledPoint> > exchangeData;
        forAll(neiProcs, i)
            exchangeData.insert
            (
                std::make_pair(neiProcs[i], LongList<labelledPoint>())
            );

        //- add points into the map
        forAll(bufferLayerPoints, pI)
        {
            const label pointI = bufferLayerPoints[pI];

            if( !(updateType[pointI] & SMOOTH) )
                continue;

            forAllRow(pProcs, pointI, i)
            {
                const label neiProc = pProcs(pointI, i);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                exchangeData[neiProc].append
                (
                    labelledPoint(globalPointLabel[pointI], pts[pointI])
                );
            }
        }

        LongList<labelledPoint> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const labelledPoint& lp = receivedData[i];

            const label triPointI = globalToLocal[lp.pointLabel()];
            pts[triPointI] = lp.coordinates();

            forAllRow(pointFacets, triPointI, ptI)
            {
                const labelledTri& tri = surf_[pointFacets(triPointI, ptI)];

                if( pointType_[tri[2]] & FACECENTRE )
                    updateType[tri[2]] |= FACECENTRE;
            }
        }
    }

    //- update coordinates of FACECENTRE vertices
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 20)
    # endif
    forAll(updateType, pI)
    {
        if( updateType[pI] & FACECENTRE )
        {
            point centre(vector::zero);
            scalar faceArea(0.0);
            forAllRow(pointFacets, pI, ptI)
            {
                const labelledTri& tri = surf_[pointFacets(pI, ptI)];
                point c(vector::zero);
                for(label i=0;i<3;++i)
                    c += pts[tri[i]];
                c /= 3;
                const scalar area = tri.mag(pts) + VSMALL;

                centre += c * area;
                faceArea += area;
            }

            pts[pI] = centre / faceArea;
        }
    }

    if( Pstream::parRun() )
        updateBufferLayers();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
