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

#include "meshOctreeAddressing.H"
#include "meshOctree.H"
#include "demandDrivenData.H"
#include "helperFunctions.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGAddressing

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeAddressing::calcGlobalPointLabels() const
{
    if( !Pstream::parRun() )
        FatalErrorIn("void meshOctreeAddressing::calcGlobalPointLabels() const")
            << "Cannot calculate global labels! Exitting" << exit(FatalError);

    const VRWGraph& nodeLabels = this->nodeLabels();
    const FRWGraph<label, 8>& nodeLeaves = this->nodeLeaves();

    //- allocate containers
    globalPointLabelPtr_ = new labelLongList(nodeLeaves.size(), -1);
    labelLongList& globalPointLabel = *globalPointLabelPtr_;

    globalPointToLocalPtr_ = new Map<label>();
    Map<label>& globalToLocal = *globalPointToLocalPtr_;

    pointProcsPtr_ = new VRWGraph(nodeLeaves.size());
    VRWGraph& pointProcs = *pointProcsPtr_;

    //- find the number of points local to each processor
    //- The point is taken to be local if it belongs to one processor, only,
    //- or to the leaf with the smallest label
    labelList nLocalPoints(Pstream::nProcs(), 0);

    forAll(nodeLeaves, pointI)
    {
        DynList<label> procs;
        procs.append(Pstream::myProcNo());
        FixedList<bool, 8> validLeaf(true);

        for(label nlI=0;nlI<8;++nlI)
        {
            const label leafI = nodeLeaves(pointI, nlI);

            if( leafI < 0 )
            {
                validLeaf[nlI] = false;
                continue;
            }

            for(label i=0;i<nlI;++i)
                if( nodeLeaves(pointI, i) == nodeLeaves(pointI, nlI) )
                {
                    validLeaf[nlI] = false;
                    validLeaf[i] = false;
                }

            procs.appendIfNotIn(octree_.returnLeaf(leafI).procNo());
        }

        label minLeaf(octree_.numberOfLeaves());
        bool found(false);

        for(label nlI=0;nlI<8;++nlI)
        {
            if( !validLeaf[nlI] )
                continue;

            const label leafI = nodeLeaves(pointI, nlI);

            minLeaf = Foam::min(leafI, minLeaf);
            found = true;
        }

        if( found && octree_.returnLeaf(minLeaf).procNo() == Pstream::myProcNo() )
        {
            if( procs.size() > 1 )
                pointProcs.setRow(pointI, procs);

            globalPointLabel[pointI] = -2;
            ++nLocalPoints[Pstream::myProcNo()];
        }
    }

    //- exchange data with other processors
    Pstream::gatherList(nLocalPoints);
    Pstream::scatterList(nLocalPoints);

    //- find the starting point label
    label startPoint(0);
    for(label procI=0;procI<Pstream::myProcNo();++procI)
        startPoint += nLocalPoints[procI];

    //- assign global labels to local points
    forAll(globalPointLabel, pointI)
    {
        if( globalPointLabel[pointI] == -2 )
        {
            globalPointLabel[pointI] = startPoint++;

            if( pointProcs.sizeOfRow(pointI) != 0 )
                globalToLocal.insert(globalPointLabel[pointI], pointI);
        }
    }

    //- distribute the labels to other processors
    //- it is done by sending the global leaf label and the node labels
    //- to processors which contain the leaves as part of buffer layers
    //- it is performed in reduce-like manner
    const labelLongList& globalLeafLabel = this->globalLeafLabel();
    const Map<label>& globalToLocalLeaf = this->globalToLocalLeafAddressing();
    const VRWGraph& leafAtProcs = this->leafAtProcs();
    const labelList& neiProcs = octree_.neiProcs();

    DynList<label> below, above;
    forAll(neiProcs, i)
    {
        if( neiProcs[i] < Pstream::myProcNo() )
        {
            above.append(neiProcs[i]);
        }
        else if( neiProcs[i] > Pstream::myProcNo() )
        {
            below.append(neiProcs[i]);
        }
    }

    VRWGraph procLeaves;
    procLeaves.reverseAddressing(Pstream::nProcs(), leafAtProcs);

    //- scatter the data from the processors above to the processors below
    //- receive the data from the processors above
    forAll(above, aboveI)
    {
        //- receive the data
        labelList receivedLabels;
        IPstream fromOtherProc(Pstream::commsTypes::blocking, above[aboveI]);
        fromOtherProc >> receivedLabels;

        label counter(0);
        while( counter < receivedLabels.size() )
        {
            const label leafI = globalToLocalLeaf[receivedLabels[counter++]];

            if( nodeLabels.sizeOfRow(leafI) == 0 )
                FatalErrorIn
                (
                   "void meshOctreeAddressing::"
                    "calcGlobalPointLabels() const"
                ) << "1. Leaf " << leafI << " is not used in the mesh!"
                 << " Exitting.." << abort(FatalError);

            for(label i=0;i<8;++i)
            {
                const label nI = nodeLabels(leafI, i);

                const label globalLabel = receivedLabels[counter++];

                const label nProcs = receivedLabels[counter++];
                for(label ppI=0;ppI<nProcs;++ppI)
                    pointProcs.appendIfNotIn(nI, receivedLabels[counter++]);

                if( globalLabel < 0 )
                    continue;

                label& gpl = globalPointLabel[nI];

                if( (gpl != -1) && (gpl != globalLabel) )
                    FatalErrorIn
                    (
                        "void meshOctreeAddressing::"
                        "calcGlobalPointLabels() const"
                    ) << "Wrong global label for point " << nI
                      << " Exitting.." << abort(FatalError);

                gpl = globalLabel;
                globalToLocal.insert(globalLabel, nI);
            }
        }
    }

    //- send the data to the processors below
    forAll(below, belowI)
    {
        const label neiProc = below[belowI];

        labelLongList dts;
        forAllRow(procLeaves, neiProc, i)
        {
            const label leafI = procLeaves(neiProc, i);

            if( nodeLabels.sizeOfRow(leafI) == 0 )
                continue;

            dts.append(globalLeafLabel[leafI]);
            for(label nI=0;nI<8;++nI)
            {
                const label nodeI = nodeLabels(leafI, nI);
                dts.append(globalPointLabel[nodeI]);

                //- add the current processor
                pointProcs.appendIfNotIn(nodeI, Pstream::myProcNo());

                dts.append(pointProcs.sizeOfRow(nodeI));
                forAllRow(pointProcs, nodeI, ppI)
                    dts.append(pointProcs(nodeI, ppI));
            }
        }

        //- send the data
        OPstream toOtherProc(Pstream::commsTypes::blocking, neiProc, dts.byteSize());
        toOtherProc << dts;
    }

    //- gather the data from the processors below to the processors above
    //- receive the data from the processors below
    forAllReverse(below, belowI)
    {
        //- receive the data
        labelList receivedLabels;
        IPstream fromOtherProc(Pstream::commsTypes::blocking, below[belowI]);
        fromOtherProc >> receivedLabels;

        label counter(0);
        while( counter < receivedLabels.size() )
        {
            const label leafI = globalToLocalLeaf[receivedLabels[counter++]];

            if( nodeLabels.sizeOfRow(leafI) == 0 )
                FatalErrorIn
                (
                   "void meshOctreeAddressing::"
                    "calcGlobalPointLabels() const"
                ) << "2. Leaf " << leafI << " is not used in the mesh!"
                 << " Exitting.." << abort(FatalError);

            for(label i=0;i<8;++i)
            {
                const label nI = nodeLabels(leafI, i);

                const label globalLabel = receivedLabels[counter++];

                const label nProcs = receivedLabels[counter++];
                for(label ppI=0;ppI<nProcs;++ppI)
                    pointProcs.appendIfNotIn(nI, receivedLabels[counter++]);

                if( globalLabel < 0 )
                    continue;

                label & gpl = globalPointLabel[nI];

                if( (gpl != -1) && (gpl != globalLabel) )
                    FatalErrorIn
                    (
                        "void meshOctreeAddressing::"
                        "calcGlobalPointLabels() const"
                    ) << "Wrong global label for point " << nI
                      << " Exitting.." << abort(FatalError);

                gpl = globalLabel;
                globalToLocal.insert(globalLabel, nI);
            }
        }
    }

    //- send the data to the processors below
    forAllReverse(above, aboveI)
    {
        const label neiProc = above[aboveI];

        labelLongList dts;
        forAllRow(procLeaves, neiProc, i)
        {
            const label leafI = procLeaves(neiProc, i);

            if( nodeLabels.sizeOfRow(leafI) == 0 )
                continue;

            dts.append(globalLeafLabel[leafI]);
            for(label nI=0;nI<8;++nI)
            {
                const label nodeI = nodeLabels(leafI, nI);
                dts.append(globalPointLabel[nodeI]);

                //- add the current processor
                pointProcs.appendIfNotIn(nodeI, Pstream::myProcNo());

                dts.append(pointProcs.sizeOfRow(nodeI));
                forAllRow(pointProcs, nodeI, ppI)
                    dts.append(pointProcs(nodeI, ppI));
            }
        }

        //- send the data
        OPstream toOtherProc(Pstream::commsTypes::blocking, neiProc, dts.byteSize());
        toOtherProc << dts;
    }
}

void meshOctreeAddressing::calcGlobalFaceLabels() const
{
    if( !Pstream::parRun() )
        FatalErrorIn("void meshOctreeAddressing::calcGlobalFaceLabels() const")
            << "Cannot calculate global labels! Exitting" << exit(FatalError);

    FatalError << "Not implemented" << exit(FatalError);
}

void meshOctreeAddressing::calcGlobalLeafLabels() const
{
    if( !Pstream::parRun() )
        FatalErrorIn("void meshOctreeAddressing::calcGlobalLeafLabels() const")
            << "Cannot calculate global labels! Exitting" << exit(FatalError);

    //- allocate the memory
    globalLeafLabelPtr_ = new labelLongList(octree_.numberOfLeaves(), -1);
    labelLongList& globalLeafLabel = *globalLeafLabelPtr_;

    globalLeafToLocalPtr_ = new Map<label>();
    Map<label>& globalToLocal = *globalLeafToLocalPtr_;

    leafAtProcsPtr_ = new VRWGraph(octree_.numberOfLeaves());
    VRWGraph& leafAtProcs = *leafAtProcsPtr_;

    //- find the number of leaves local to each processor
    labelList nLeavesAtProc(Pstream::nProcs(), 0);

    label nLeaves(0);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(static) reduction(+:nLeaves)
    # endif
    for(label leafI=0;leafI<octree_.numberOfLeaves();++leafI)
    {
        const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);

        if( oc.procNo() == Pstream::myProcNo() )
            ++nLeaves;
    }

    nLeavesAtProc[Pstream::myProcNo()] = nLeaves;

    //- exchange the data with other processors
    Pstream::gatherList(nLeavesAtProc);
    Pstream::scatterList(nLeavesAtProc);

    //- find the starting labels for
    nLeaves = 0;
    for(label procI=0;procI<Pstream::myProcNo();++procI)
        nLeaves += nLeavesAtProc[procI];

    //- set the global labels to local leaves
    labelLongList otherProcLeaves;
    for(label leafI=0;leafI<octree_.numberOfLeaves();++leafI)
    {
        const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);

        if( oc.procNo() == Pstream::myProcNo() )
        {
            globalLeafLabel[leafI] = nLeaves++;
        }
        else
        {
            otherProcLeaves.append(leafI);
            leafAtProcs.append(leafI, Pstream::myProcNo());
            leafAtProcs.append(leafI, oc.procNo());
        }
    }

    //- the rest of the code is needed in case an additional layer of
    //- of octree leaves belonging to other processors is added in order to
    //- simplify the procedure for generation of mesh templates

    //- allocate the map for exchanging of data
    std::map<label, LongList<meshOctreeCubeBasic> > exchangeData;
    const labelList& neiProcs = octree_.neiProcs();
    forAll(neiProcs, procI)
    {
        const std::pair<label, LongList<meshOctreeCubeBasic> > pp
        (
            neiProcs[procI],
            LongList<meshOctreeCubeBasic>()
        );

        exchangeData.insert(pp);
    }

    //- Here we have to combine the information from all processors
    //- it is started such that all processors send the leaves to the processor
    //- that contains them locally

    //- fill the map with data
    forAll(otherProcLeaves, i)
    {
        const meshOctreeCubeBasic& oc = octree_.returnLeaf(otherProcLeaves[i]);
        meshOctreeCubeBasic coc(oc);
        coc.setProcNo(Pstream::myProcNo());
        exchangeData[oc.procNo()].append(coc);
    }

    //- exchange the data with other processors
    LongList<meshOctreeCubeBasic> rLeaves;
    help::exchangeMap(exchangeData, rLeaves, Pstream::commsTypes::scheduled);

    //- update the local data
    forAll(rLeaves, i)
    {
        const label cLabel = octree_.findLeafLabelForPosition(rLeaves[i]);

        globalToLocal.insert(globalLeafLabel[cLabel], cLabel);
        leafAtProcs.appendIfNotIn(cLabel, Pstream::myProcNo());
        leafAtProcs.appendIfNotIn(cLabel, rLeaves[i].procNo());
    }

    //- now the global leaf labels shall be sent from the processors
    //- that own the leaves to the processors that also contain them
    std::map<label, labelLongList> exchangeLabels;
    std::map<label, LongList<meshOctreeCubeBasic> >::iterator it;
    for(it=exchangeData.begin();it!=exchangeData.end();++it)
    {
        it->second.clear();
        exchangeLabels.insert(std::make_pair(it->first, labelLongList()));
    }

    //- fill in the data
    forAll(leafAtProcs, leafI)
    {
        if( octree_.returnLeaf(leafI).procNo() == Pstream::myProcNo() )
        {
            forAllRow(leafAtProcs, leafI, i)
            {
                const label procI = leafAtProcs(leafI, i);

                if( procI == Pstream::myProcNo() )
                    continue;

                exchangeData[procI].append(octree_.returnLeaf(leafI));
                exchangeLabels[procI].append(globalLeafLabel[leafI]);
            }
        }
    }

    //- exchange the data
    rLeaves.clear();
    help::exchangeMap(exchangeData, rLeaves, Pstream::commsTypes::scheduled);
    labelLongList rLabels;
    help::exchangeMap(exchangeLabels, rLabels, Pstream::commsTypes::scheduled);

    if( rLeaves.size() != rLabels.size() )
        FatalErrorIn("void meshOctreeAddressing::calcGlobalLeafLabels() const")
            << "Invalid list size!" << abort(FatalError);

    //- set the labels to the leaves originating from other processors
    forAll(rLeaves, i)
    {
        const label cLabel = octree_.findLeafLabelForPosition(rLeaves[i]);

        globalLeafLabel[cLabel] = rLabels[i];
        globalToLocal.insert(rLabels[i], cLabel);
    }

    //- update leafAtProcs for all processors
    exchangeLabels.clear();
    forAll(neiProcs, procI)
        exchangeLabels.insert(std::make_pair(neiProcs[procI], labelLongList()));

    forAllConstIter(Map<label>, globalToLocal, iter)
    {
        const label leafI = iter();

        if( octree_.returnLeaf(leafI).procNo() != Pstream::myProcNo() )
            continue;

        forAllRow(leafAtProcs, leafI, i)
        {
            const label procI = leafAtProcs(leafI, i);

            if( procI == Pstream::myProcNo() )
                continue;

            labelLongList& dts = exchangeLabels[procI];
            dts.append(iter.key());

            dts.append(leafAtProcs.sizeOfRow(leafI));
            forAllRow(leafAtProcs, leafI, j)
                dts.append(leafAtProcs(leafI, j));
        }
    }

    //- exchange the data
    rLabels.clear();
    help::exchangeMap(exchangeLabels, rLabels, Pstream::commsTypes::scheduled);

    //- update the local data
    label counter(0);
    while( counter < rLabels.size() )
    {
        const label gLabel = rLabels[counter++];

        if( !globalToLocal.found(gLabel) )
            FatalError << "Cannot find global label " << gLabel
                << exit(FatalError);

        const label leafI = globalToLocal[gLabel];

        const label numberOfProcs = rLabels[counter++];
        for(label i=0;i<numberOfProcs;++i)
            leafAtProcs.append(leafI, rLabels[counter++]);
    }

    # ifdef DEBUGAddressing
    returnReduce(1, sumOp<label>());
    const List<direction>& boxType = this->boxType();
    const VRWGraph& nodeLabels = this->nodeLabels();
    for(label i=0;i<Pstream::nProcs();++i)
    {
        if( i == Pstream::myProcNo() )
        {
            forAll(globalLeafLabel, leafI)
            {
                Pout << "Leaf " << leafI << "of type " << label(boxType[leafI])
                    << " has global label " << globalLeafLabel[leafI]
                    << " and coordinates " << octree_.returnLeaf(leafI)
                    << " and located at procs " << leafAtProcs[leafI]
                    << " node labels " << nodeLabels[leafI] << endl;

                if( octree_.returnLeaf(leafI).procNo() == Pstream::myProcNo() )
                    continue;
                if( globalToLocal[globalLeafLabel[leafI]] != leafI )
                    FatalError << "Crap!!" << abort(FatalError);
            }
        }

        returnReduce(1, sumOp<label>());
    }

    //- check if the leaf at procs is ok
    exchangeData.clear();
    for(label i=0;i<Pstream::nProcs();++i)
    {
        if( i == Pstream::myProcNo() )
            continue;

        exchangeData[i];
    }

    for(label leafI=0;leafI<octree_.numberOfLeaves();++leafI)
        for(label i=0;i<Pstream::nProcs();++i)
        {
            if( i == Pstream::myProcNo() )
                continue;

            exchangeData[i].append(octree_.returnLeaf(leafI));
        }

    std::map<label, List<meshOctreeCubeBasic> > rMap;
    help::exchangeMap(exchangeData, rMap);

    for(std::map<label, List<meshOctreeCubeBasic> >::const_iterator it=rMap.begin();it!=rMap.end();++it)
    {
        const List<meshOctreeCubeBasic>& data = it->second;

        forAll(data, i)
        {
            const label leafI = octree_.findLeafLabelForPosition(data[i]);

            if( leafI < 0 )
                continue;

            if( !leafAtProcs.contains(leafI, it->first) )
                FatalError << "Problem!!" << leafI
                    << " does not contain processor " << it->first
                    << " contains procs " << leafAtProcs[leafI]
                    << abort(FatalError);
        }
    }

    returnReduce(1, sumOp<label>());
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
