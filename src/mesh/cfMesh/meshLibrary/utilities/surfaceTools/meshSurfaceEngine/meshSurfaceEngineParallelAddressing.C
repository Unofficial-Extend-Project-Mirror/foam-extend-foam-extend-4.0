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

#include "meshSurfaceEngine.H"
#include "demandDrivenData.H"
#include "Map.H"
#include "HashSet.H"
#include "helperFunctionsPar.H"

#include <map>
#include <set>

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceEngine::calcGlobalBoundaryPointLabels() const
{
    if( !globalBoundaryPointLabelPtr_ )
        globalBoundaryPointLabelPtr_ = new labelList();

    const labelList& bPoints = this->boundaryPoints();
    labelList& globalPointLabel = *globalBoundaryPointLabelPtr_;
    globalPointLabel.setSize(bPoints.size());
    globalPointLabel = -1;

    if( !bpProcsPtr_ )
        bpProcsPtr_ = new VRWGraph(bPoints.size());

    if( !globalBoundaryPointToLocalPtr_ )
        globalBoundaryPointToLocalPtr_ = new Map<label>();

    if( !bpNeiProcsPtr_ )
        bpNeiProcsPtr_ = new DynList<label>();

    if( !Pstream::parRun() )
        return;

    VRWGraph& bpAtProcs = *bpProcsPtr_;

    //- find local points for the given processor
    const faceListPMG& faces = mesh_.faces();
    const labelList& bp = this->bp();
    const PtrList<processorBoundaryPatch>& procBoundaries =
        mesh_.procBoundaries();

    //- find which processor contain a given bnd point
    forAll(procBoundaries, patchI)
    {
        const label start = procBoundaries[patchI].patchStart();
        const label end = start + procBoundaries[patchI].patchSize();
        for(label faceI=start;faceI<end;++faceI)
        {
            const face& f = faces[faceI];
            forAll(f, pI)
                if( bp[f[pI]] != -1 )
                {
                    bpAtProcs.appendIfNotIn
                    (
                        bp[f[pI]],
                        procBoundaries[patchI].myProcNo()
                    );
                    bpAtProcs.appendIfNotIn
                    (
                        bp[f[pI]],
                        procBoundaries[patchI].neiProcNo()
                    );
                }
        }
    }

    //- exchange data with other processor and update bpAtProcs
    //- continue this process as long as there is some change
    bool finished;
    do
    {
        finished = true;

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label end = start + procBoundaries[patchI].patchSize();

            labelLongList dts;
            labelHashSet addedPoint;
            for(label faceI=start;faceI<end;++faceI)
            {
                const face& f = faces[faceI];
                forAll(f, pI)
                    if( (bp[f[pI]] != -1) && !addedPoint.found(f[pI]) )
                    {
                        addedPoint.insert(f[pI]);
                        const label bpI = bp[f[pI]];
                        //- data is sent as follows
                        //- 1. face position in patch
                        //- 2. local point position in face
                        //- 3. number of processors for point
                        //- 4. proc labels
                        dts.append(faceI-start);
                        dts.append((f.size()-pI)%f.size());
                        dts.append(bpAtProcs.sizeOfRow(bpI));
                        forAllRow(bpAtProcs, bpI, i)
                            dts.append(bpAtProcs(bpI, i));
                    }
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                dts.byteSize()
            );
            toOtherProc << dts;
        }

        forAll(procBoundaries, patchI)
        {
            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo()
            );
            labelList receivedData;
            fromOtherProc >> receivedData;

            const label start = procBoundaries[patchI].patchStart();

            label counter(0);
            while( counter < receivedData.size() )
            {
                const face& f = faces[start+receivedData[counter++]];
                const label pI = receivedData[counter++];
                const label nProcs = receivedData[counter++];
                for(label i=0;i<nProcs;++i)
                {
                    const label neiProc = receivedData[counter++];
                    if( !bpAtProcs.contains(bp[f[pI]], neiProc) )
                    {
                        bpAtProcs.append(bp[f[pI]], neiProc);
                        finished = false;
                    }
                }
            }
        }

        reduce(finished, minOp<bool>());
    } while( !finished );

    //- start calculating global point labels
    label nLocalPoints(0);
    boolList localPoints(bPoints.size(), true);
    forAll(bpAtProcs, bpI)
        if( bpAtProcs.sizeOfRow(bpI) != 0 )
        {
            label pMin = bpAtProcs(bpI, 0);
            forAllRow(bpAtProcs, bpI, procI)
                if( bpAtProcs(bpI, procI) < pMin )
                    pMin = bpAtProcs(bpI, procI);

            if( pMin == Pstream::myProcNo() )
            {
                ++nLocalPoints;
            }
            else
            {
                localPoints[bpI] = false;
            }
        }
        else
        {
            ++nLocalPoints;
        }

    //- give local points their labels
    label startPoint(0);
    labelList nPointsAtProc(Pstream::nProcs());
    nPointsAtProc[Pstream::myProcNo()] = nLocalPoints;
    Pstream::gatherList(nPointsAtProc);
    Pstream::scatterList(nPointsAtProc);
    for(label i=0;i<Pstream::myProcNo();++i)
        startPoint += nPointsAtProc[i];

    forAll(localPoints, pI)
        if( localPoints[pI] )
            globalPointLabel[pI] = startPoint++;

    //- send labels to non-local points
    do
    {
        finished = true;

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label end = start + procBoundaries[patchI].patchSize();

            labelLongList dts;
            labelHashSet addedPoint;
            for(label faceI=start;faceI<end;++faceI)
            {
                const face& f = faces[faceI];
                forAll(f, pI)
                {
                    const label bpI = bp[f[pI]];
                    if( (bpI != -1) && (globalPointLabel[bpI] != -1) )
                    {
                        if( addedPoint.found(bpI) )
                            continue;

                        addedPoint.insert(bpI);
                        //- data is sent as follows
                        //- 1. face position in patch
                        //- 2. local point position in face
                        //- 3. global point label
                        dts.append(faceI-start);
                        dts.append((f.size()-pI)%f.size());
                        dts.append(globalPointLabel[bpI]);
                    }
                }
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                dts.byteSize()
            );
            toOtherProc << dts;
        }

        forAll(procBoundaries, patchI)
        {
            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo()
            );
            labelList receivedData;
            fromOtherProc >> receivedData;

            const label start = procBoundaries[patchI].patchStart();

            label counter(0);
            while( counter < receivedData.size() )
            {
                const face& f = faces[start+receivedData[counter++]];
                const label pI = receivedData[counter++];
                const label globalLabel = receivedData[counter++];

                if( globalPointLabel[bp[f[pI]]] == -1 )
                {
                    globalPointLabel[bp[f[pI]]] = globalLabel;
                    finished = false;
                }
                else if( globalPointLabel[bp[f[pI]]] != globalLabel )
                {
                    FatalErrorIn
                    (
                        "void meshSurfaceEngine::"
                        "calcGlobalBoundaryPointLabels() const"
                    ) << "Point labels in proc boundary "
                        << procBoundaries[patchI].patchName()
                        << " face " << f << " pI = " << pI
                        << nl << " label " << globalPointLabel[bp[f[pI]]]
                        << nl << " other global label " << globalLabel
                        << " do not match!" << abort(FatalError);
                }
            }
        }

        reduce(finished, minOp<bool>());
    } while( !finished );

    //- create globalToLocal addressing
    forAll(bpAtProcs, bpI)
    {
        if( bpAtProcs.sizeOfRow(bpI) != 0 )
        {
            globalBoundaryPointToLocalPtr_->insert(globalPointLabel[bpI], bpI);

            forAllRow(bpAtProcs, bpI, j)
            {
                const label procI = bpAtProcs(bpI, j);

                if( procI == Pstream::myProcNo() )
                    continue;

                bpNeiProcsPtr_->appendIfNotIn(procI);
            }
        }
    }
}

void meshSurfaceEngine::calcGlobalBoundaryEdgeLabels() const
{
    if( !globalBoundaryEdgeLabelPtr_ )
        globalBoundaryEdgeLabelPtr_ = new labelList();

    const edgeList& bEdges = this->edges();

    labelList& globalEdgeLabel = *globalBoundaryEdgeLabelPtr_;
    globalEdgeLabel.setSize(bEdges.size());
    globalEdgeLabel = -1;

    if( !beProcsPtr_ )
        beProcsPtr_ = new VRWGraph(bEdges.size());

    if( !globalBoundaryEdgeToLocalPtr_ )
        globalBoundaryEdgeToLocalPtr_ = new Map<label>();

    if( !beNeiProcsPtr_ )
        beNeiProcsPtr_ = new DynList<label>();

    if( !Pstream::parRun() )
        return;

    VRWGraph& beAtProcs = *beProcsPtr_;

    const faceListPMG& faces = mesh_.faces();
    const labelList& bp = this->bp();
    const VRWGraph& pEdges = this->boundaryPointEdges();
    const PtrList<processorBoundaryPatch>& procBoundaries =
        mesh_.procBoundaries();

    //- find which processors contain a given bnd edge
    typedef std::set<std::pair<label, label> > procEdgeMap;
    List<procEdgeMap> facesWithProcBndEdges(procBoundaries.size());

    forAll(procBoundaries, patchI)
    {
        const label start = procBoundaries[patchI].patchStart();
        const label end = start + procBoundaries[patchI].patchSize();
        for(label faceI=start;faceI<end;++faceI)
        {
            const face& f = faces[faceI];

            forAll(f, eI)
            {
                if( bp[f[eI]] != -1 )
                {
                    const edge e = f.faceEdge(eI);

                    const label bpI = bp[e.start()];
                    label edgeI(-1);
                    forAllRow(pEdges, bpI, peI)
                        if( bEdges[pEdges(bpI, peI)] == e )
                        {
                            edgeI = pEdges(bpI, peI);
                            break;
                        }

                    if( edgeI != -1 )
                    {
                        facesWithProcBndEdges[patchI].insert
                        (
                            std::make_pair(faceI, eI)
                        );

                        beAtProcs.appendIfNotIn
                        (
                            edgeI,
                            procBoundaries[patchI].myProcNo()
                        );
                        beAtProcs.appendIfNotIn
                        (
                            edgeI,
                            procBoundaries[patchI].neiProcNo()
                        );
                    }
                }
            }
        }
    }

    //- exchange data with other processor and update beAtProcs
    //- continue this process as long as there is some change
    bool finished;
    do
    {
        finished = true;

        forAll(facesWithProcBndEdges, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const procEdgeMap& procBndEdges = facesWithProcBndEdges[patchI];

            labelLongList dts;
            forAllConstIter(procEdgeMap, procBndEdges, it)
            {
                const std::pair<label, label>& fPair = *it;

                const face& f = faces[fPair.first];
                const edge e = f.faceEdge(fPair.second);

                if( bp[e.start()] != -1 )
                {
                    const label bpI = bp[e.start()];
                    label edgeI(-1);
                    forAllRow(pEdges, bpI, peI)
                        if( bEdges[pEdges(bpI, peI)] == e )
                        {
                            edgeI = pEdges(bpI, peI);
                            break;
                        }

                    if( edgeI != -1 )
                    {
                        //- data is sent as follows
                        //- 1. face position in patch
                        //- 2. local edge position in face
                        //- 3. number of processors for edge
                        //- 4. proc labels
                        dts.append(fPair.first-start);
                        dts.append((f.size()-fPair.second-1)%f.size());
                        dts.append(beAtProcs.sizeOfRow(edgeI));
                        forAllRow(beAtProcs, edgeI, i)
                            dts.append(beAtProcs(edgeI, i));
                    }
                }
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                dts.byteSize()
            );
            toOtherProc << dts;
        }

        forAll(procBoundaries, patchI)
        {
            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo()
            );
            labelList receivedData;
            fromOtherProc >> receivedData;

            const label start = procBoundaries[patchI].patchStart();

            label counter(0);
            while( counter < receivedData.size() )
            {
                const label faceI = start+receivedData[counter++];
                const face& f = faces[faceI];
                const label eI = receivedData[counter++];

                const edge e = f.faceEdge(eI);
                label edgeI(-1);
                forAllRow(pEdges, bp[e.start()], peI)
                    if( bEdges[pEdges(bp[e.start()], peI)] == e )
                        edgeI = pEdges(bp[e.start()], peI);

                const label nProcs = receivedData[counter++];
                for(label i=0;i<nProcs;++i)
                {
                    const label neiProc = receivedData[counter++];
                    if( !beAtProcs.contains(edgeI, neiProc) )
                    {
                        facesWithProcBndEdges[patchI].insert
                        (
                            std::make_pair(faceI, eI)
                        );

                        beAtProcs.append(edgeI, neiProc);
                        finished = false;
                    }
                }
            }
        }

        reduce(finished, minOp<bool>());
    } while( !finished );

    //- start calculating global edge labels
    label nLocalEdges(0);
    boolList localEdges(bEdges.size(), true);
    forAll(beAtProcs, bpI)
        if( beAtProcs.sizeOfRow(bpI) != 0 )
        {
            label pMin = beAtProcs(bpI, 0);
            forAllRow(beAtProcs, bpI, procI)
                if( beAtProcs(bpI, procI) < pMin )
                    pMin = beAtProcs(bpI, procI);

            if( pMin == Pstream::myProcNo() )
            {
                ++nLocalEdges;
            }
            else
            {
                localEdges[bpI] = false;
            }
        }
        else
        {
            ++nLocalEdges;
        }

    //- give local points their labels
    label startEdge(0);
    labelList nEdgesAtProc(Pstream::nProcs());
    nEdgesAtProc[Pstream::myProcNo()] = nLocalEdges;
    Pstream::gatherList(nEdgesAtProc);
    Pstream::scatterList(nEdgesAtProc);
    for(label i=0;i<Pstream::myProcNo();++i)
        startEdge += nEdgesAtProc[i];

    forAll(localEdges, pI)
        if( localEdges[pI] )
            globalEdgeLabel[pI] = startEdge++;

    //- send labels to non-local edges
    do
    {
        finished = true;

        forAll(facesWithProcBndEdges, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const procEdgeMap& procBndEdges = facesWithProcBndEdges[patchI];

            labelLongList dts;
            forAllConstIter(procEdgeMap, procBndEdges, it)
            {
                const label faceI = it->first;
                const face& f = faces[faceI];

                const label eI = it->second;
                const edge e = f.faceEdge(eI);

                if( bp[e.start()] != -1 )
                {
                    const label bpI = bp[e.start()];
                    label edgeI(-1);
                    forAllRow(pEdges, bpI, peI)
                        if( bEdges[pEdges(bpI, peI)] == e )
                        {
                            edgeI = pEdges(bpI, peI);
                            break;
                        }

                    if( (edgeI != -1) && (globalEdgeLabel[edgeI] != -1) )
                    {
                        //- data is sent as follows
                        //- 1. face position in patch
                        //- 2. local edge position in face
                        //- 3. global edge label
                        dts.append(faceI-start);
                        dts.append((f.size()-eI-1)%f.size());
                        dts.append(globalEdgeLabel[edgeI]);
                    }
                }
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                dts.byteSize()
            );
            toOtherProc << dts;
        }

        forAll(procBoundaries, patchI)
        {
            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo()
            );
            labelList receivedData;
            fromOtherProc >> receivedData;

            const label start = procBoundaries[patchI].patchStart();

            label counter(0);
            while( counter < receivedData.size() )
            {
                const face& f = faces[start+receivedData[counter++]];
                const label eI = receivedData[counter++];

                const edge e = f.faceEdge(eI);
                label edgeI(-1);
                forAllRow(pEdges, bp[e.start()], peI)
                    if( bEdges[pEdges(bp[e.start()], peI)] == e )
                        edgeI = pEdges(bp[e.start()], peI);

                const label globalLabel = receivedData[counter++];

                if( globalEdgeLabel[edgeI] == -1 )
                {
                    globalEdgeLabel[edgeI] = globalLabel;
                    finished = false;
                }
                else if( globalEdgeLabel[edgeI] != globalLabel )
                {
                    FatalErrorIn
                    (
                        "void meshSurfaceEngine::"
                        "calcGlobalBoundaryEdgeLabels() const"
                    ) << "Edge labels do not match!" << abort(FatalError);
                }
            }
        }

        reduce(finished, minOp<bool>());
    } while( !finished );

    //- create globalToLocal addressing
    forAll(beAtProcs, beI)
    {
        if( beAtProcs.sizeOfRow(beI) != 0 )
        {
            globalBoundaryEdgeToLocalPtr_->insert(globalEdgeLabel[beI], beI);

            forAllRow(beAtProcs, beI, j)
            {
                const label procI = beAtProcs(beI, j);

                if( procI == Pstream::myProcNo() )
                    continue;

                beNeiProcsPtr_->appendIfNotIn(procI);
            }
        }
    }
}

void meshSurfaceEngine::calcAddressingForProcEdges() const
{
    const labelList& globalEdgeLabel = this->globalBoundaryEdgeLabel();
    const labelList& boundaryFacePatches = this->boundaryFacePatches();
    const VRWGraph& eFaces = this->edgeFaces();
    const VRWGraph& beAtProcs = this->beAtProcs();
    const Map<label>& globalToLocal = this->globalToLocalBndEdgeAddressing();
    const DynList<label>& beNeiProcs = this->beNeiProcs();

    std::map<label, labelLongList> exchangeData;
    forAll(beNeiProcs, i)
        exchangeData.insert(std::make_pair(beNeiProcs[i], labelLongList()));

    //- check if it the surface is manifold over inter-processor edges
    Map<label> nFacesAtEdge;
    forAllConstIter(Map<label>, globalToLocal, iter)
    {
        const label beI = iter();
        nFacesAtEdge.insert(beI, eFaces.sizeOfRow(beI));

        forAllRow(beAtProcs, beI, i)
        {
            const label neiProc = beAtProcs(beI, i);

            if( neiProc == Pstream::myProcNo() )
                continue;

            labelLongList& dts = exchangeData[neiProc];
            dts.append(iter.key());
            dts.append(eFaces.sizeOfRow(beI));
        }
    }

    labelLongList receivedData;
    help::exchangeMap(exchangeData, receivedData);
    for(label counter=0;counter<receivedData.size();)
    {
        const label beI = globalToLocal[receivedData[counter++]];
        nFacesAtEdge[beI] += receivedData[counter++];
    }

    forAllConstIter(Map<label>, nFacesAtEdge, iter)
    {
        if( iter() != 2 )
            FatalErrorIn
            (
                "void meshSurfaceEngine::calcAddressingForProcEdges() const"
            ) << "Surface is not manifold at boundary edge "
              << iter.key() << exit(FatalError);
    }

    //- find the processor which contains other face containing an edge
    //- at inter-processor boundary
    exchangeData.clear();
    forAll(beNeiProcs, i)
        exchangeData.insert(std::make_pair(beNeiProcs[i], labelLongList()));

    forAllConstIter(Map<label>, globalToLocal, iter)
    {
        const label beI = iter();

        forAllRow(beAtProcs, beI, procI)
        {
            const label neiProc = beAtProcs(beI, procI);

            if( neiProc == Pstream::myProcNo() )
                continue;

            if( eFaces.sizeOfRow(beI) == 1 )
            {
                exchangeData[neiProc].append(globalEdgeLabel[beI]);
                exchangeData[neiProc].append
                (
                    boundaryFacePatches[eFaces(beI, 0)]
                );
            }
        }
    }

    //- exchange edge-face patches with other processors
    std::map<label, labelList> receivedMap;
    help::exchangeMap(exchangeData, receivedMap);

    //- store other face patches in a map
    otherEdgeFaceAtProcPtr_ = new Map<label>();
    otherEdgeFacePatchPtr_ = new Map<label>();
    Map<label>& otherProcPatches = *otherEdgeFacePatchPtr_;
    Map<label>& otherFaceProc = *otherEdgeFaceAtProcPtr_;
    for
    (
        std::map<label, labelList>::const_iterator iter=receivedMap.begin();
        iter!=receivedMap.end();
        ++iter
    )
    {
        const labelList& receivedData = iter->second;

        label counter(0);
        while( counter < receivedData.size() )
        {
            const label beI = globalToLocal[receivedData[counter++]];
            const label patch = receivedData[counter++];
            if( eFaces.sizeOfRow(beI) == 1 )
            {
                otherProcPatches.insert(beI, patch);
                otherFaceProc.insert(beI, iter->first);
            }
        }
    }
}

void meshSurfaceEngine::calcGlobalBoundaryFaceLabels() const
{
    const faceList::subList& bFaces = boundaryFaces();

    if( !globalBoundaryFaceLabelPtr_ )
        globalBoundaryFaceLabelPtr_ = new labelList(bFaces.size());

    labelList& globalFaceLabel = *globalBoundaryFaceLabelPtr_;

    labelList nFacesAtProc(Pstream::nProcs());
    nFacesAtProc[Pstream::myProcNo()] = bFaces.size();
    Pstream::gatherList(nFacesAtProc);
    Pstream::scatterList(nFacesAtProc);

    label startFace(0);
    for(label i=0;i<Pstream::myProcNo();++i)
        startFace += nFacesAtProc[i];

    forAll(bFaces, fI)
        globalFaceLabel[fI] = startFace++;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
