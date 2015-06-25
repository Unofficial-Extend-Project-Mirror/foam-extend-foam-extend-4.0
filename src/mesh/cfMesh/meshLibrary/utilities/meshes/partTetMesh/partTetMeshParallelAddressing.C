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
#include "polyMeshGenModifier.H"
#include "partTetMesh.H"
#include "polyMeshGenAddressing.H"
#include "helperFunctionsPar.H"
#include "parPartTet.H"

#include <map>

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void partTetMesh::createParallelAddressing
(
    const labelLongList& nodeLabelForPoint,
    const labelLongList& nodeLabelForFace,
    const labelLongList& nodeLabelForCell
)
{
    //- vertices marked as SMOOTH and BOUNDARY are used by the smoother
    const direction useType = SMOOTH + BOUNDARY;
    
    //- allocate global point labels
    if( !globalPointLabelPtr_ )
        globalPointLabelPtr_ = new labelLongList();
    labelLongList& globalTetPointLabel = *globalPointLabelPtr_;
    globalTetPointLabel.setSize(points_.size());
    globalTetPointLabel = -1;
    
    //- allocated point-processors addressing
    if( !pAtProcsPtr_ )
        pAtProcsPtr_ = new VRWGraph();
    VRWGraph& pProcs = *pAtProcsPtr_;
    pProcs.setSize(0);
    pProcs.setSize(points_.size());
    
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
    
    const polyMeshGenAddressing& addressing = origMesh_.addressingData();
    const Map<label>& globalToLocalPointAddressing =
        addressing.globalToLocalPointAddressing();
    const VRWGraph& pAtProcs = addressing.pointAtProcs();
    const DynList<label>& pNeiProcs = addressing.pointNeiProcs();
    
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
            !smoothVertex_[nodeLabelForPoint[pI]]
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
        
        smoothVertex_[nodeLabelForPoint[pointI]] = NONE;
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
        if( !(smoothVertex_[nodeLabelForPoint[pI]] & useType) )
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
    nPointsAtProc[Pstream::myProcNo()] = points_.size() - nSharedPoints;
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
        
        if( !(smoothVertex_[pLabel] & useType) )
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
        
        globalTetPointLabel[pLabel] = startPoint++;
        
        forAllRow(pAtProcs, pI, procI)
        {
            const label neiProc = pAtProcs(pI, procI);
            
            if( neiProc == Pstream::myProcNo() )
                continue;
            
            //- the following information is sent to other processor
            //- 1. global point label in the original mesh
            //- 2. global point label in the tet mesh
            exchangeData[neiProc].append(it.key());
            exchangeData[neiProc].append(globalTetPointLabel[pLabel]);
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
        
        globalTetPointLabel[pLabel] = tgI;
    }
    
    //- set global labels for remaining points
    forAll(globalTetPointLabel, pI)
    {
        if( globalTetPointLabel[pI] == -1 )
            globalTetPointLabel[pI] = startPoint++;
    }
        
    //- create global to local mapping
    forAll(globalTetPointLabel, pI)
    {
        if( pProcs.sizeOfRow(pI) != 0 )
        {
            pAtParallelBoundaries.append(pI);
            globalToLocal.insert(globalTetPointLabel[pI], pI);
        }
    }
    
    //- mark vertices at parallel boundaries
    forAll(smoothVertex_, pI)
        if( (smoothVertex_[pI] & useType) && (pProcs.sizeOfRow(pI) != 0) )
            smoothVertex_[pI] |= PARALLELBOUNDARY;
        
    //- create neighbour processors addressing
    if( !neiProcsPtr_ )
        neiProcsPtr_ = new DynList<label>();
    DynList<label>& neiProcs = *neiProcsPtr_;
    
    for(iter=exchangeData.begin();iter!=exchangeData.end();++iter)
        neiProcs.append(iter->first);
    
    # ifdef DEBUGSmooth
    for(label i=0;i<Pstream::nProcs();++i)
    {
        if( i == Pstream::myProcNo() )
        {
            Pout << "globalTetPointLabel " << globalTetPointLabel << endl;
        }
        
        returnReduce(i, sumOp<label>());
    }
    
    returnReduce(1, sumOp<label>());
    
    forAll(nodeLabelForPoint, pI)
    {
        const label tpI = nodeLabelForPoint[pI];
        if( tpI != -1 && globalTetPointLabel[tpI] == -1 )
            FatalError << "Crap1 " << tpI << abort(FatalError);
    }
    
    returnReduce(1, sumOp<label>());
    
    forAll(nodeLabelForFace, fI)
    {
        const label tpI = nodeLabelForFace[fI];
        if( tpI != -1 && globalTetPointLabel[tpI] == -1 )
        {
            Pout << "Face point " << tpI << " is at procs "
                << pProcs[tpI] << endl;
            FatalError << "Crap2" << tpI << abort(FatalError);
        }
    }
    
    returnReduce(1, sumOp<label>());
    
    forAll(nodeLabelForCell, cI)
    {
        const label tpI = nodeLabelForCell[cI];
        if( tpI != -1 && globalTetPointLabel[tpI] == -1 )
            FatalError << "Crap3" << tpI << abort(FatalError);
    }
    
    forAll(smoothVertex_, vI)
        if( smoothVertex_[vI] & partTetMesh::PARALLELBOUNDARY )
            Pout << "Point " << globalTetPointLabel[vI]
            << " is at par bnd" << endl;
        
    Serr << Pstream::myProcNo() << "points " << points_ << endl;
    Serr << Pstream::myProcNo() << "Tets " << tets_ << endl;
    forAll(pProcs, pI)
    {
        if( pProcs.sizeOfRow(pI) == 0 )
            continue;
        
        Serr << Pstream::myProcNo() << "Point " << globalTetPointLabel[pI]
            << " is at procs " << pProcs[pI] << " n tets "
            << pointTets_[pI].size() << endl;
    }
    
    returnReduce(1, sumOp<label>());
    # endif
}

void partTetMesh::createBufferLayers()
{
    VRWGraph& pProcs = *pAtProcsPtr_;
    labelLongList& globalTetPointLabel = *globalPointLabelPtr_;
    Map<label>& globalToLocal = *globalToLocalPointAddressingPtr_;
    const DynList<label>& neiProcs = *this->neiProcsPtr_;
    
    if( !pAtBufferLayersPtr_ )
        pAtBufferLayersPtr_ = new labelLongList();
    labelLongList& pAtBufferLayers = *pAtBufferLayersPtr_;
    pAtBufferLayers.clear();
    
    //- create the map
    std::map<label, LongList<parPartTet> > exchangeTets;
    forAll(neiProcs, procI)
        exchangeTets.insert
        (
            std::make_pair(neiProcs[procI], LongList<parPartTet>())
        );
    
    //- go through the tets and add the ones having vertices at parallel
    //- boundaries for sending
    forAll(tets_, tetI)
    {
        const partTet& pt = tets_[tetI];
        
        DynList<label> sendToProcs;
        forAll(pt, i)
        {
            const label pLabel = pt[i];
            
            if( smoothVertex_[pLabel] & PARALLELBOUNDARY )
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
            const parPartTet tet
            (
                labelledPoint(globalTetPointLabel[pt[0]], points_[pt[0]]),
                labelledPoint(globalTetPointLabel[pt[1]], points_[pt[1]]),
                labelledPoint(globalTetPointLabel[pt[2]], points_[pt[2]]),
                labelledPoint(globalTetPointLabel[pt[3]], points_[pt[3]])
            );
            
            forAll(sendToProcs, i)
            {
                exchangeTets[sendToProcs[i]].append(tet);
                
                forAll(pt, j)
                {
                    if( pProcs.sizeOfRow(pt[j]) == 0 )
                        pAtBufferLayers.append(pt[j]);
                    
                    pProcs.appendIfNotIn(pt[j], sendToProcs[i]);
                }
            }
        }
    }
    
    LongList<parPartTet> receivedTets;
    help::exchangeMap(exchangeTets, receivedTets);
    exchangeTets.clear();
    
    Map<label> newGlobalToLocal;
    forAll(receivedTets, i)
    {
        const parPartTet& tet = receivedTets[i];
        
        DynList<label> tetPointLabels;
        for(label j=0;j<4;++j)
        {
            const label gpI = tet[j].pointLabel();
            
            if( globalToLocal.found(gpI) )
            {
                const label pI = globalToLocal[gpI];
                pointTets_.append(pI, tets_.size());
                tetPointLabels.append(pI);
            }
            else if( newGlobalToLocal.found(gpI) )
            {
                tetPointLabels.append(newGlobalToLocal[gpI]);
            }
            else
            {
                newGlobalToLocal.insert(gpI, points_.size());
                tetPointLabels.append(points_.size());
                points_.append(tet[j].coordinates());
                nodeLabelInOrigMesh_.append(-1);
                smoothVertex_.append(NONE);
                DynList<label> helper;
                helper.append(tets_.size());
                pointTets_.appendList(helper);
                
                globalTetPointLabel.append(gpI);
                helper[0] = Pstream::myProcNo();
                pProcs.appendList(helper);
            }
        }
        
        //- append tet
        tets_.append
        (
            partTet
            (
                tetPointLabels[0],
                tetPointLabels[1],
                tetPointLabels[2],
                tetPointLabels[3]
            )
        );
    }
    
    //- insert the global labels of the buffer points
    //- into the globalToLocal map
    forAllConstIter(Map<label>, newGlobalToLocal, it)
        globalToLocal.insert(it.key(), it());

    # ifdef DEBUGSmooth
    for(label i=0;i<Pstream::nProcs();++i)
    {
        if( Pstream::myProcNo() == i )
        {
            forAllConstIter(Map<label>, globalToLocal, it)
            {
                DynList<label> np;
                if( it() < pProcs.size() )
                forAllRow(pProcs, it(), j)
                    np.append(pProcs(it(), j));
            
                Pout << "Tet mesh point " << it() << " has global label "
                    << it.key() << " and is located at procs "
                    << np << endl;
            }
        }
        
        returnReduce(1, sumOp<label>());
    }
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
