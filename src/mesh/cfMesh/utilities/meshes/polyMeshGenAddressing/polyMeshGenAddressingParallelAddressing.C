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

\*---------------------------------------------------------------------------*/

#include "polyMeshGenAddressing.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenAddressing::calcGlobalPointLabels() const
{
    if( !Pstream::parRun() )
        FatalErrorIn
        (
            "void polyMeshGenAddressing::calcGlobalPointLabels() const"
        ) << "Cannot calculate global point addressing "
            << "for a serial run!" << exit(FatalError);

    if( !globalPointLabelPtr_ )
        globalPointLabelPtr_ = new labelLongList();
    labelLongList& globalPointLabel = *globalPointLabelPtr_;
    globalPointLabel.setSize(mesh_.points().size());
    globalPointLabel = -1;

    if( !pProcsPtr_ )
        pProcsPtr_ = new VRWGraph();
    VRWGraph& pProcs = *pProcsPtr_;
    pProcs.setSize(mesh_.points().size());

    if( !globalToLocalPointAddressingPtr_ )
        globalToLocalPointAddressingPtr_ = new Map<label>();
    Map<label>& globalToLocal = *globalToLocalPointAddressingPtr_;

    if( !pointNeiProcsPtr_ )
        pointNeiProcsPtr_ = new DynList<label>();
    DynList<label>& pNeiProcs = *pointNeiProcsPtr_;

    const faceListPMG& faces = mesh_.faces();
    const PtrList<processorBoundaryPatch>& procBoundaries =
        mesh_.procBoundaries();

    //- find which processors contain a given bnd point
    List<std::map<label, std::pair<label, label> > > patchPoints;
    patchPoints.setSize(procBoundaries.size());
    forAll(procBoundaries, patchI)
    {
        const label start = procBoundaries[patchI].patchStart();
        const label end = start + procBoundaries[patchI].patchSize();

        std::map<label, std::pair<label, label> >& patchPointsMap =
            patchPoints[patchI];

        for(label faceI=start;faceI<end;++faceI)
        {
            const face& f = faces[faceI];

            forAll(f, pI)
            {
                std::map<label, std::pair<label, label> >::iterator it =
                    patchPointsMap.find(f[pI]);

                if( it != patchPointsMap.end() )
                    continue;

                const std::pair<label, label> pp
                (
                    faceI-start,
                    (f.size()-pI)%f.size()
                );
                patchPointsMap.insert(std::make_pair(f[pI], pp));

                pProcs.appendIfNotIn(f[pI], procBoundaries[patchI].myProcNo());
                pProcs.appendIfNotIn(f[pI], procBoundaries[patchI].neiProcNo());
            }
        }
    }

    //- exchange data with other processor and update bProcs
    //- continue this process as long as there is some change
    bool finished;
    do
    {
        finished = true;

        forAll(procBoundaries, patchI)
        {
            const std::map<label, std::pair<label, label> >& patchPointsMap =
                patchPoints[patchI];
            std::map<label, std::pair<label, label> >::const_iterator it;

            labelLongList dataToSend;
            for(it=patchPointsMap.begin();it!=patchPointsMap.end();++it)
            {
                //- data is sent as follows
                //- 1. face position in patch
                //- 2. local point position in face
                //- 3. number of patches for point
                //- 4. patch labels
                dataToSend.append(it->second.first);
                dataToSend.append(it->second.second);
                dataToSend.append(pProcs.sizeOfRow(it->first));
                forAllRow(pProcs, it->first, i)
                    dataToSend.append(pProcs(it->first, i));
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                dataToSend.byteSize()
            );
            toOtherProc << dataToSend;
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
                    if( !pProcs.contains(f[pI], neiProc) )
                    {
                        pProcs.append(f[pI], neiProc);
                        finished = false;
                    }
                }
            }
        }

        reduce(finished, minOp<bool>());
    } while( !finished );

    //- start calculating global point labels
    label nLocalPoints(0);
    boolList localPoints(mesh_.points().size(), true);
    forAll(pProcs, pI)
        if( pProcs.sizeOfRow(pI) != 0 )
        {
            label pMin = pProcs(pI, 0);
            forAllRow(pProcs, pI, procI)
                pMin = Foam::min(pMin, pProcs(pI, procI));

            if( pMin == Pstream::myProcNo() )
            {
                ++nLocalPoints;
            }
            else
            {
                localPoints[pI] = false;
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
            const std::map<label, std::pair<label, label> >& patchPointsMap =
                patchPoints[patchI];
            std::map<label, std::pair<label, label> >::const_iterator it;

            labelLongList dataToSend;
            for(it=patchPointsMap.begin();it!=patchPointsMap.end();++it)
            {
                if( globalPointLabel[it->first] != -1 )
                {
                    //- data is sent as follows
                    //- 1. face position in patch
                    //- 2. local point position in face
                    //- 3. global point label
                    dataToSend.append(it->second.first);
                    dataToSend.append(it->second.second);
                    dataToSend.append(globalPointLabel[it->first]);
                }
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                dataToSend.byteSize()
            );
            toOtherProc << dataToSend;
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

                if( globalPointLabel[f[pI]] == -1 )
                {
                    globalPointLabel[f[pI]] = globalLabel;
                    finished = false;
                }
                else if( globalPointLabel[f[pI]] != globalLabel )
                {
                    FatalErrorIn
                    (
                        "void polyMeshGenAddressing::"
                        "calcGlobalPointLabels() const"
                    ) << "Point labels on processors do not match!"
                        << abort(FatalError);
                }
            }
        }

        reduce(finished, minOp<bool>());
    } while( !finished );

    //- create global to local addressing and pointNeiProcs addressing
    forAll(globalPointLabel, pointI)
        if( pProcs.sizeOfRow(pointI) != 0 )
        {
            globalToLocal.insert(globalPointLabel[pointI], pointI);

            forAllRow(pProcs, pointI, i)
                pNeiProcs.appendIfNotIn(pProcs(pointI, i));
        }
}

void polyMeshGenAddressing::calcGlobalFaceLabels() const
{
    if( !globalFaceLabelPtr_ )
        globalFaceLabelPtr_ = new labelLongList();
    labelLongList& globalFaceLabel = *globalFaceLabelPtr_;
    globalFaceLabel.setSize(mesh_.faces().size());
    globalFaceLabel = -1;

    if( !Pstream::parRun() )
        return;

    const label nIntFaces = mesh_.nInternalFaces();

    const PtrList<processorBoundaryPatch>& procBoundaries =
        mesh_.procBoundaries();

    label nBoundaryFaces = 0;
    forAll(procBoundaries, patchI)
    {
        if( procBoundaries[patchI].owner() )
            nBoundaryFaces += procBoundaries[patchI].patchSize();
    }

    label startFace(0);
    labelList numberOfInternalFaces(Pstream::nProcs());
    numberOfInternalFaces[Pstream::myProcNo()] = nIntFaces + nBoundaryFaces;

    Pstream::gatherList(numberOfInternalFaces);
    Pstream::scatterList(numberOfInternalFaces);
    for(label i=0;i<Pstream::myProcNo();++i)
        startFace += numberOfInternalFaces[i];

    //- calculate labels for internal faces
    for(label faceI=0;faceI<nIntFaces;++faceI)
        globalFaceLabel[faceI] = startFace++;

    //- calculate labels for processor boundaries
    forAll(procBoundaries, patchI)
    {
        if( procBoundaries[patchI].owner() )
        {
            const label start = procBoundaries[patchI].patchStart();
            labelList dataToSend(procBoundaries[patchI].patchSize());
            forAll(dataToSend, bfI)
            {
                globalFaceLabel[start+bfI] = startFace;
                dataToSend[bfI] = startFace++;
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                dataToSend.byteSize()
            );
            toOtherProc << dataToSend;
        }
    }

    forAll(procBoundaries, patchI)
    {
        if( !procBoundaries[patchI].owner() )
        {
            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            labelList receivedLabels;
            fromOtherProc >> receivedLabels;

            const label start = procBoundaries[patchI].patchStart();
            forAll(receivedLabels, bfI)
            {
                globalFaceLabel[start+bfI] = receivedLabels[bfI];
            }
        }
    }

    //- create global labels for boundary faces
    startFace = 0;
    forAll(numberOfInternalFaces, procI)
        startFace += numberOfInternalFaces[procI];

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    labelListList numberOfBoundaryFaces(Pstream::nProcs());
    numberOfBoundaryFaces[Pstream::myProcNo()].setSize(boundaries.size());
    forAll(boundaries, patchI)
    {
        numberOfBoundaryFaces[Pstream::myProcNo()][patchI] =
            boundaries[patchI].patchSize();
    }
    Pstream::gatherList(numberOfBoundaryFaces);
    Pstream::scatterList(numberOfBoundaryFaces);

    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label nBoundaryFaces = boundaries[patchI].patchSize();

        for(label procI=0;procI<Pstream::myProcNo();++procI)
            startFace += numberOfBoundaryFaces[procI][patchI];

        for(label bfI=0;bfI<nBoundaryFaces;++bfI)
            globalFaceLabel[start+bfI] = startFace++;

        for(label procI=Pstream::myProcNo()+1;procI<Pstream::nProcs();++procI)
            startFace += numberOfBoundaryFaces[procI][patchI];
    }
}

void polyMeshGenAddressing::calcGlobalCellLabels() const
{
    if( !globalCellLabelPtr_ )
        globalCellLabelPtr_ = new labelLongList();
    labelLongList& globalCellLabel = *globalCellLabelPtr_;
    globalCellLabel.setSize(mesh_.cells().size());
    globalCellLabel = -1;

    if( !Pstream::parRun() )
        return;

    label startLabel(0);
    labelList nCellsAtProc(Pstream::nProcs());
    nCellsAtProc[Pstream::myProcNo()] = globalCellLabel.size();
    Pstream::gatherList(nCellsAtProc);
    Pstream::scatterList(nCellsAtProc);
    for(label i=0;i<Pstream::myProcNo();++i)
        startLabel += nCellsAtProc[i];

    forAll(globalCellLabel, cellI)
        globalCellLabel[cellI] = startLabel++;
}

void polyMeshGenAddressing::calcGlobalEdgeLabels() const
{
    if( !Pstream::parRun() )
        FatalErrorIn
        (
            "void polyMeshGenAddressing::calcGlobalEdgeLabels() const"
        ) << "Cannot calculate global edge addressing "
            << "for a serial run?!?!" << exit(FatalError);

    if( !globalEdgeLabelPtr_ )
        globalEdgeLabelPtr_ = new labelLongList();
    labelLongList& globalEdgeLabel = *globalEdgeLabelPtr_;
    globalEdgeLabel.setSize(this->edges().size());
    globalEdgeLabel = -1;

    if( !eProcsPtr_ )
        eProcsPtr_ = new VRWGraph();
    VRWGraph& eProcs = *eProcsPtr_;
    eProcs.setSize(this->edges().size());

    if( !globalToLocalEdgeAddressingPtr_ )
        globalToLocalEdgeAddressingPtr_ = new Map<label>();
    Map<label>& globalToLocal = *globalToLocalEdgeAddressingPtr_;

    if( !edgeNeiProcsPtr_ )
        edgeNeiProcsPtr_ = new DynList<label>();
    DynList<label>& eNeiProcs = *edgeNeiProcsPtr_;

    const faceListPMG& faces = mesh_.faces();
    const VRWGraph& faceEdges = this->faceEdges();
    const PtrList<processorBoundaryPatch>& procBoundaries =
        mesh_.procBoundaries();

    //- find which processor contain a given bnd edge
    forAll(procBoundaries, patchI)
    {
        const label start = procBoundaries[patchI].patchStart();
        const label end = start + procBoundaries[patchI].patchSize();
        for(label faceI=start;faceI<end;++faceI)
        {
            forAllRow(faceEdges, faceI, eI)
            {
                eProcs.appendIfNotIn
                (
                    faceEdges(faceI, eI),
                    procBoundaries[patchI].myProcNo()
                );
                eProcs.appendIfNotIn
                (
                    faceEdges(faceI, eI),
                    procBoundaries[patchI].neiProcNo()
                );
            }
        }
    }

    //- exchange data with other processor and update eProcs
    //- continue this process as long as there is some change
    bool finished;
    do
    {
        finished = true;

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label end = start + procBoundaries[patchI].patchSize();

/*            label nToSend(0);
            for(label faceI=start;faceI<end;++faceI)
            {
                forAllRow(faceEdges, faceI, eI)
                {
                    nToSend += 3;
                    nToSend += eProcs.sizeOfRow(faceEdges(faceI, eI));
                }
            }
*/
            labelLongList dataToSend;
            //nToSend = 0;
            for(label faceI=start;faceI<end;++faceI)
            {
                forAllRow(faceEdges, faceI, eI)
                {
                    const label edgeI = faceEdges(faceI, eI);
                    const face& f = faces[faceI];
                    //- data is sent as follows
                    //- 1. face position in patch
                    //- 2. local edge position in face
                    //- 3. number of processors for edge
                    //- 4. processor labels
                    dataToSend.append(faceI-start);
                    dataToSend.append((f.size()-eI-1)%f.size());
                    dataToSend.append(eProcs.sizeOfRow(edgeI));
                    forAllRow(eProcs, edgeI, i)
                        dataToSend.append(eProcs(edgeI, i));
                }
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                dataToSend.byteSize()
            );
            toOtherProc << dataToSend;
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
                const label eI = receivedData[counter++];

                const label edgeI = faceEdges(faceI, eI);

                const label nProcs = receivedData[counter++];
                for(label i=0;i<nProcs;++i)
                {
                    const label neiProc = receivedData[counter++];
                    if( !eProcs.contains(edgeI, neiProc) )
                    {
                        eProcs.append(edgeI, neiProc);
                        finished = false;
                    }
                }
            }
        }

        reduce(finished, minOp<bool>());
    } while( !finished );

    //- start calculating global edge labels
    label nLocalEdges(0);
    boolList localEdges(eProcs.size(), true);
    forAll(eProcs, edgeI)
        if( eProcs.sizeOfRow(edgeI) != 0 )
        {
            label pMin = eProcs(edgeI, 0);
            forAllRow(eProcs, edgeI, procI)
                if( eProcs(edgeI, procI) < pMin )
                    pMin = eProcs(edgeI, procI);

            if( pMin == Pstream::myProcNo() )
            {
                ++nLocalEdges;
            }
            else
            {
                localEdges[edgeI] = false;
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

        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label end = start + procBoundaries[patchI].patchSize();

            labelLongList dataToSend;
            for(label faceI=start;faceI<end;++faceI)
            {
                forAllRow(faceEdges, faceI, eI)
                {
                    const label edgeI = faceEdges(faceI, eI);
                    if( globalEdgeLabel[edgeI] != -1 )
                    {
                        const face& f = faces[faceI];
                        //- data is sent as follows
                        //- 1. face position in patch
                        //- 2. local edge position in face
                        //- 3. number global label for the edge
                        dataToSend.append(faceI-start);
                        dataToSend.append((f.size()-eI-1)%f.size());
                        dataToSend.append(globalEdgeLabel[edgeI]);
                    }
                }
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                dataToSend.byteSize()
            );
            toOtherProc << dataToSend;
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
                const label eI = receivedData[counter++];

                const label edgeI = faceEdges(faceI, eI);

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
                        "void polyMeshGenAddressing::"
                        "calcGlobalEdgeLabels() const"
                    ) << "Edge labels on processors do not match!"
                        << abort(FatalError);
                }
            }
        }

        reduce(finished, minOp<bool>());
    } while( !finished );

    //- create globalToLocal addressing
    forAll(eProcs, edgeI)
        if( eProcs.sizeOfRow(edgeI) != 0 )
        {
            globalToLocal.insert(globalEdgeLabel[edgeI], edgeI);

            forAllRow(eProcs, edgeI, i)
                eNeiProcs.appendIfNotIn(eProcs(edgeI, i));
        }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelLongList& polyMeshGenAddressing::globalPointLabel() const
{
    if( !globalPointLabelPtr_ || !pProcsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const labelLongList& polyMeshGenAddressing::"
                "globalPointLabel() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcGlobalPointLabels();
    }

    return *globalPointLabelPtr_;
}

const labelLongList& polyMeshGenAddressing::globalFaceLabel() const
{
    if( !globalFaceLabelPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const labelLongList& polyMeshGenAddressing"
                "::globalFaceLabel() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcGlobalFaceLabels();
    }

    return *globalFaceLabelPtr_;
}

const labelLongList& polyMeshGenAddressing::globalCellLabel() const
{
    if( !globalCellLabelPtr_ )
        calcGlobalCellLabels();

    return *globalCellLabelPtr_;
}

const labelLongList& polyMeshGenAddressing::globalEdgeLabel() const
{
    if( !globalEdgeLabelPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const labelLongList& polyMeshGenAddressing"
                "::globalEdgeLabel() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcGlobalEdgeLabels();
    }

    return *globalEdgeLabelPtr_;
}

const VRWGraph& polyMeshGenAddressing::pointAtProcs() const
{
    if( !globalPointLabelPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& polyMeshGenAddressing::pointAtProcs() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcGlobalPointLabels();
    }

    return *pProcsPtr_;
}

const DynList<label>& polyMeshGenAddressing::pointNeiProcs() const
{
    if( !pointNeiProcsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const DynList<label>& polyMeshGenAddressing"
                "::pointNeiProcs() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcGlobalPointLabels();
    }

    return *pointNeiProcsPtr_;
}

const Map<label>& polyMeshGenAddressing::globalToLocalPointAddressing() const
{
    if( !globalToLocalPointAddressingPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const Map<label>& polyMeshGenAddressing"
                "::globalToLocalPointAddressing() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcGlobalPointLabels();
    }

    return *globalToLocalPointAddressingPtr_;
}

const VRWGraph& polyMeshGenAddressing::edgeAtProcs() const
{
    if( !globalEdgeLabelPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& polyMeshGenAddressing::edgeAtProcs() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcGlobalEdgeLabels();
    }

    return *eProcsPtr_;
}

const DynList<label>& polyMeshGenAddressing::edgeNeiProcs() const
{
    if( !edgeNeiProcsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const DynList<label>& polyMeshGenAddressing::edgeNeiProcs() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcGlobalEdgeLabels();
    }

    return *edgeNeiProcsPtr_;
}

const Map<label>& polyMeshGenAddressing::globalToLocalEdgeAddressing() const
{
    if( !globalToLocalEdgeAddressingPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const Map<label>& polyMeshGenAddressing"
                "::globalToLocalEdgeAddressing() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcGlobalEdgeLabels();
    }

    return *globalToLocalEdgeAddressingPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
