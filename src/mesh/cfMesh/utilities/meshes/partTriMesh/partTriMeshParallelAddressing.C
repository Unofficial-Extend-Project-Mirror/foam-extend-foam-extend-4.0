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
#include "meshSurfacePartitioner.H"
#include "partTriMesh.H"
#include "triSurfModifier.H"
#include "meshSurfaceEngine.H"
#include "helperFunctionsPar.H"
#include "parTriFace.H"

#include <map>

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void partTriMesh::createParallelAddressing
(
    const labelList& nodeLabelForPoint,
    const labelList& nodeLabelForFace
)
{
    const meshSurfaceEngine& mse = mPart_.surfaceEngine();

    const pointField& pts = surf_.points();

    //- vertices marked as SMOOTH are used by the smoother
    const direction useType = SMOOTH;

    //- allocate global point labels
    if( !globalPointLabelPtr_ )
        globalPointLabelPtr_ = new labelLongList();
    labelLongList& globalPointLabel = *globalPointLabelPtr_;
    globalPointLabel.setSize(pts.size());
    globalPointLabel = -1;

    //- allocated point-processors addressing
    if( !pAtProcsPtr_ )
        pAtProcsPtr_ = new VRWGraph();
    VRWGraph& pProcs = *pAtProcsPtr_;
    pProcs.setSize(0);
    pProcs.setSize(pts.size());

    //- allocate global-to-local point addressing
    if( !globalToLocalPointAddressingPtr_ )
        globalToLocalPointAddressingPtr_ = new Map<label>();
    Map<label>& globalToLocal = *globalToLocalPointAddressingPtr_;
    globalToLocal.clear();

    //- allocate storage for points at parallel boundaries
    if( !pAtParallelBoundariesPtr_ )
        pAtParallelBoundariesPtr_ = new labelLongList();
    labelLongList& pAtParallelBoundaries = *pAtParallelBoundariesPtr_;
    pAtParallelBoundaries.clear();

    //- create point-processors addressing
    std::map<label, labelLongList> exchangeData;
    std::map<label, labelLongList>::iterator iter;

    const Map<label>& globalToLocalPointAddressing =
        mse.globalToLocalBndPointAddressing();
    const VRWGraph& pAtProcs = mse.bpAtProcs();
    const DynList<label>& pNeiProcs = mse.bpNeiProcs();

    forAll(pNeiProcs, procI)
        exchangeData.insert(std::make_pair(pNeiProcs[procI], labelLongList()));

    //- make sure that the same vertices are marked for smoothing on all procs
    //- this is performed by sending the labels of vertices which are not used
    //- for tet mesh creation and the tet mesh vertices which are not moved
    forAllConstIter(Map<label>, globalToLocalPointAddressing, it)
    {
        const label pI = it();

        if(
            nodeLabelForPoint[pI] == -1 ||
            !pointType_[nodeLabelForPoint[pI]]
        )
        {
            forAllRow(pAtProcs, pI, procI)
            {
                const label neiProc = pAtProcs(pI, procI);
                if( neiProc == Pstream::myProcNo() )
                    continue;

                exchangeData[neiProc].append(it.key());
            }
        }
    }

    //- exchange data with other processors
    labelLongList receivedData;
    help::exchangeMap(exchangeData, receivedData);

    //- set the values according to other processors
    forAll(receivedData, i)
    {
        const label pointI = globalToLocalPointAddressing[receivedData[i]];

        if( nodeLabelForPoint[pointI] == -1 )
            continue;

        pointType_[nodeLabelForPoint[pointI]] = NONE;
    }

    for(iter=exchangeData.begin();iter!=exchangeData.end();++iter)
        iter->second.clear();

    //- start creating global-to-local addressing
    //- find the starting point labels
    label startPoint(0), nLocalPoints(0), nSharedPoints(0);

    //- count the number of points at processor boundaries
    forAllConstIter(Map<label>, globalToLocalPointAddressing, it)
    {
        const label pI = it();

        if( nodeLabelForPoint[pI] == -1 )
            continue;
        if( !(pointType_[nodeLabelForPoint[pI]] & useType) )
            continue;

        ++nSharedPoints;

        label pMin(Pstream::myProcNo());
        forAllRow(pAtProcs, pI, procI)
            pMin = Foam::min(pMin, pAtProcs(pI, procI));

        if( pMin == Pstream::myProcNo() )
            ++nLocalPoints;
    }

    labelList nPointsAtProc(Pstream::nProcs());
    nSharedPoints -= nLocalPoints;
    nPointsAtProc[Pstream::myProcNo()] = pts.size() - nSharedPoints;
    Pstream::gatherList(nPointsAtProc);
    Pstream::scatterList(nPointsAtProc);

    for(label i=0;i<Pstream::myProcNo();++i)
        startPoint += nPointsAtProc[i];

    //- create global labels for points at processor boundaries
    forAllConstIter(Map<label>, globalToLocalPointAddressing, it)
    {
        const label pI = it();

        if( nodeLabelForPoint[pI] == -1 )
            continue;

        const label pLabel = nodeLabelForPoint[pI];

        if( !(pointType_[pLabel] & useType) )
            continue;

        label pMin(Pstream::myProcNo());
        forAllRow(pAtProcs, pI, procI)
        {
            const label neiProc = pAtProcs(pI, procI);
            pProcs.append(pLabel, neiProc);
            pMin = Foam::min(pMin, neiProc);
        }

        if( pMin != Pstream::myProcNo() )
            continue;

        globalPointLabel[pLabel] = startPoint++;

        forAllRow(pAtProcs, pI, procI)
        {
            const label neiProc = pAtProcs(pI, procI);

            if( neiProc == Pstream::myProcNo() )
                continue;

            //- the following information is sent to other processor
            //- 1. global point label in the original mesh
            //- 2. global point label in the tet mesh
            exchangeData[neiProc].append(it.key());
            exchangeData[neiProc].append(globalPointLabel[pLabel]);
        }
    }

    //- exchange data with other processors
    receivedData.clear();
    help::exchangeMap(exchangeData, receivedData);

    label counter(0);
    while( counter < receivedData.size() )
    {
        const label gpI = receivedData[counter++];
        const label tgI = receivedData[counter++];
        const label pLabel =
            nodeLabelForPoint[globalToLocalPointAddressing[gpI]];

        globalPointLabel[pLabel] = tgI;
    }

    //- set global labels for remaining points
    forAll(globalPointLabel, pI)
    {
        if( globalPointLabel[pI] == -1 )
            globalPointLabel[pI] = startPoint++;
    }

    //- create global to local mapping
    forAll(globalPointLabel, pI)
    {
        if( pProcs.sizeOfRow(pI) != 0 )
        {
            pAtParallelBoundaries.append(pI);
            globalToLocal.insert(globalPointLabel[pI], pI);
        }
    }

    //- mark vertices at parallel boundaries
    forAll(pointType_, pI)
        if( (pointType_[pI] & useType) && (pProcs.sizeOfRow(pI) != 0) )
            pointType_[pI] |= PARALLELBOUNDARY;

    //- create neighbour processors addressing
    if( !neiProcsPtr_ )
        neiProcsPtr_ = new DynList<label>();
    DynList<label>& neiProcs = *neiProcsPtr_;

    for(iter=exchangeData.begin();iter!=exchangeData.end();++iter)
        neiProcs.append(iter->first);
}

void partTriMesh::createBufferLayers()
{
    pointField& pts = triSurfModifier(surf_).pointsAccess();

    VRWGraph& pProcs = *pAtProcsPtr_;
    labelLongList& globalPointLabel = *globalPointLabelPtr_;
    Map<label>& globalToLocal = *globalToLocalPointAddressingPtr_;
    const DynList<label>& neiProcs = *this->neiProcsPtr_;

    if( !pAtBufferLayersPtr_ )
        pAtBufferLayersPtr_ = new labelLongList();
    labelLongList& pAtBufferLayers = *pAtBufferLayersPtr_;
    pAtBufferLayers.clear();

    //- create the map
    std::map<label, LongList<parTriFace> > exchangeTrias;
    forAll(neiProcs, procI)
        exchangeTrias.insert
        (
            std::make_pair(neiProcs[procI], LongList<parTriFace>())
        );

    //- loop over triangles and add the ones having vertices at parallel
    //- boundaries for sending
    forAll(surf_, triI)
    {
        const labelledTri& pt = surf_[triI];

        DynList<label> sendToProcs;
        forAll(pt, i)
        {
            const label pLabel = pt[i];

            if( pointType_[pLabel] & PARALLELBOUNDARY )
            {
                forAllRow(pProcs, pLabel, i)
                {
                    const label neiProc = pProcs(pLabel, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    sendToProcs.appendIfNotIn(neiProc);
                }
            }
        }

        if( sendToProcs.size() )
        {
            const parTriFace tri
            (
                globalPointLabel[pt[0]],
                globalPointLabel[pt[1]],
                globalPointLabel[pt[2]],
                triangle<point, point>(pts[pt[0]], pts[pt[1]], pts[pt[2]])
            );

            forAll(sendToProcs, i)
            {
                exchangeTrias[sendToProcs[i]].append(tri);

                forAll(pt, j)
                {
                    if( pProcs.sizeOfRow(pt[j]) == 0 )
                        pAtBufferLayers.append(pt[j]);

                    pProcs.appendIfNotIn(pt[j], sendToProcs[i]);
                }
            }
        }
    }

    //- receive triangles sent to this processor
    std::map<label, List<parTriFace> > receivedTriangles;
    help::exchangeMap(exchangeTrias, receivedTriangles);
    exchangeTrias.clear();

    //- add triangles into the mesh and update the addressing
    Map<label> newGlobalToLocal;
    std::map<label, point> addCoordinates;
    label nPoints = pts.size();
    for
    (
        std::map<label, List<parTriFace> >::const_iterator it =
        receivedTriangles.begin();
        it!=receivedTriangles.end();
        ++it
    )
    {
        const List<parTriFace>& receivedTrias = it->second;

        forAll(receivedTrias, i)
        {
            const parTriFace& tri = receivedTrias[i];

            DynList<label, 3> triPointLabels(3);
            for(label j=0;j<3;++j)
            {
                const label gpI = tri.globalLabelOfPoint(j);

                if( globalToLocal.found(gpI) )
                {
                    //- point already exists in the triangulation
                    const label pI = globalToLocal[gpI];
                    triPointLabels[j] = pI;
                }
                else if( newGlobalToLocal.found(gpI) )
                {
                    //- point is already added into the triangulation
                    triPointLabels[j] = newGlobalToLocal[gpI];
                    pProcs.appendIfNotIn(newGlobalToLocal[gpI], it->first);
                }
                else
                {
                    //- point does not exist in the triangulation
                    //- and is not yet added in
                    newGlobalToLocal.insert(gpI, nPoints);
                    triPointLabels[j] = nPoints;

                    point tp;
                    if( j == 0 )
                    {
                        tp = tri.trianglePoints().a();
                    }
                    else if( j == 1 )
                    {
                        tp = tri.trianglePoints().b();
                    }
                    else
                    {
                        tp = tri.trianglePoints().c();
                    }
                    addCoordinates[nPoints] = tp;
                    ++nPoints;

                    pointLabelInMeshSurface_.append(-1);
                    pointType_.append(NONE);

                    DynList<label> triAtProcs;
                    triAtProcs.append(it->first);

                    globalPointLabel.append(gpI);
                    triAtProcs.append(Pstream::myProcNo());
                    pProcs.appendList(triAtProcs);
                }
            }

            //- append tet
            surf_.appendTriangle
            (
                labelledTri
                (
                    triPointLabels[0],
                    triPointLabels[1],
                    triPointLabels[2],
                    -1
                )
            );
        }
    }

    //- store newly added points
    pts.setSize(nPoints);
    for
    (
        std::map<label, point>::const_iterator it=addCoordinates.begin();
        it!=addCoordinates.end();
        ++it
    )
        pts[it->first] = it->second;

    addCoordinates.clear();

    //- insert the global labels of the buffer points
    //- into the globalToLocal map
    forAllConstIter(Map<label>, newGlobalToLocal, it)
        globalToLocal.insert(it.key(), it());

    //- update addressing of the surface mesh
    surf_.clearAddressing();
}

void partTriMesh::updateBufferLayers()
{
    const pointField& points = surf_.points();
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

        forAllRow(pProcs, pointI, i)
        {
            const label neiProc = pProcs(pointI, i);

            if( neiProc == Pstream::myProcNo() )
                continue;

            exchangeData[neiProc].append
            (
                labelledPoint(globalPointLabel[pointI], points[pointI])
            );
        }
    }

    LongList<labelledPoint> receivedData;
    help::exchangeMap(exchangeData, receivedData);

    forAll(receivedData, i)
    {
        const labelledPoint& lp = receivedData[i];

        this->updateVertex
        (
            globalToLocal[lp.pointLabel()],
            lp.coordinates()
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
