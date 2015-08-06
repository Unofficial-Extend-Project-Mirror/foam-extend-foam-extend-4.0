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

#include "meshOctreeModifier.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
# endif

//#define OCTREETiming
//#define DEBUGBalancing

# ifdef DEBUGBalancing
#include <sstream>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Private member functions

void meshOctreeModifier::loadDistribution(const direction usedType)
{
    if( octree_.neiProcs().size() == 0 )
        return;

    # ifdef OCTREETiming
    returnReduce(1, sumOp<label>());
    const scalar startTime = omp_get_wtime();

    returnReduce(1, sumOp<label>());
    const scalar t1 = omp_get_wtime();
    Info << "Creation of list of leaves lasted " << t1-startTime << endl;
    # endif

    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    label localNumWeighs(0);
    labelList globalLeafWeight(leaves.size());

    if( usedType )
    {
        forAll(leaves, leafI)
        {
            if( leaves[leafI]->cubeType() & usedType )
            {
                globalLeafWeight[leafI] = localNumWeighs;
                ++localNumWeighs;
            }
            else
            {
                globalLeafWeight[leafI] = -1;
            }
        }
    }
    else
    {
        forAll(leaves, leafI)
        {
            globalLeafWeight[leafI] = localNumWeighs;
            ++localNumWeighs;
        }
    }

    # ifdef OCTREETiming
    returnReduce(1, sumOp<label>());
    const scalar t2 = omp_get_wtime();
    Info << "Creation of global leaf weights lasted " << t2-t1 << endl;
    # endif

    const label totalNumWeights = returnReduce(localNumWeighs, sumOp<label>());
    const label nWeightsPerProcessor = totalNumWeights / Pstream::nProcs();

    //- check if balancing should be performed
    //- the tolerance is set to 5% difference in the number of boxes
    //- from the ideal one
    label doBalancing(0);
    if(
        mag
        (
            scalar(localNumWeighs - nWeightsPerProcessor) /
            nWeightsPerProcessor
        ) > 0.05
    )
        doBalancing = 1;

    reduce(doBalancing, maxOp<label>());

    if( doBalancing == 0 )
        return;

    Info << "Distributing load between processors" << endl;

    //- start calculating new partitions
    //- find global labels of the leaf boxes
    doBalancing = 0;

    labelList procWeights(Pstream::nProcs());
    procWeights[Pstream::myProcNo()] = localNumWeighs;
    Pstream::gatherList(procWeights);
    Pstream::scatterList(procWeights);

    for(label procI=0;procI<Pstream::myProcNo();++procI)
        doBalancing += procWeights[procI];

    forAll(globalLeafWeight, lI)
    {
        if( globalLeafWeight[lI] != -1 )
            globalLeafWeight[lI] += doBalancing;
    }

    //- leaf boxes which are not in the range for the current processor
    //- shall be migrated to other processors
    std::map<label, labelLongList> leavesToSend;

    bool oneRemainingBox(false);
    forAll(globalLeafWeight, leafI)
    {
        if( globalLeafWeight[leafI] == -1 )
            continue;
        if( !oneRemainingBox && (leafI == leaves.size() -1) )
            continue;

        const label newProc =
            Foam::min
            (
                globalLeafWeight[leafI] / nWeightsPerProcessor,
                Pstream::nProcs()-1
            );

        if( newProc != Pstream::myProcNo() )
        {
            leavesToSend[newProc].append(leafI);
            leaves[leafI]->setProcNo(newProc);

            # ifdef DEBUGBalancing
            if( leaves[leafI]->hasContainedElements() )
                Serr << Pstream::myProcNo() << "Deleting a DATA cube "
                << leaves[leafI]->coordinates() << " data is "
                << leaves[leafI]->containedElements() << endl;
            # endif
        }
        else
        {
            oneRemainingBox = true;
        }
    }

    # ifdef OCTREETiming
    returnReduce(1, sumOp<label>());
    const scalar t3 = omp_get_wtime();
    Info << "Completed assignment of leaves to processors in " << t3-t2 << endl;
    # endif

    //- send the information to other processors
    //- all processors shall received a list containing the same information
    //- each processor informs which other processors shall receive data from
    //- that processor
    labelListList sendToProcesssors(Pstream::nProcs());
    sendToProcesssors[Pstream::myProcNo()].setSize(leavesToSend.size());
    label counter(0);
    for
    (
        std::map<label, labelLongList>::const_iterator it=leavesToSend.begin();
        it!=leavesToSend.end();
        ++it
    )
        sendToProcesssors[Pstream::myProcNo()][counter++] = it->first;

    Pstream::gatherList(sendToProcesssors);
    Pstream::scatterList(sendToProcesssors);

    labelHashSet receiveFrom;
    forAll(sendToProcesssors, procI)
        forAll(sendToProcesssors[procI], neiI)
            if( sendToProcesssors[procI][neiI] == Pstream::myProcNo() )
                receiveFrom.insert(procI);

    //- receive coordinates from processors with lower labels
    LongList<meshOctreeCubeBasic> migratedCubes;
    forAllConstIter(labelHashSet, receiveFrom, iter)
    {
        if( iter.key() >= Pstream::myProcNo() )
            continue;

        List<meshOctreeCubeBasic> mc;

        IPstream fromOtherProc(Pstream::blocking, iter.key());

        fromOtherProc >> mc;

        label currSize = migratedCubes.size();
        migratedCubes.setSize(currSize+mc.size());
        forAll(mc, mcI)
        {
            migratedCubes[currSize] = mc[mcI];
            ++currSize;
        }
    }

    //- send the coordinates of the boxes to processors with greater label
    const labelList& sendToProcs = sendToProcesssors[Pstream::myProcNo()];
    forAll(sendToProcs, i)
    {
        const label procI = sendToProcs[i];

        if( procI <= Pstream::myProcNo() )
            continue;

        List<meshOctreeCubeBasic> sendCoordinates
        (
            leavesToSend[procI].size()
        );

        forAll(leavesToSend[procI], lI)
        {
            const meshOctreeCube& oc = *leaves[leavesToSend[procI][lI]];
            sendCoordinates[lI] =
                meshOctreeCubeBasic
                (
                    oc.coordinates(),
                    oc.cubeType()
                );
        }

        OPstream toOtherProc
        (
            Pstream::blocking,
            procI,
            sendCoordinates.byteSize()
        );

        toOtherProc << sendCoordinates;
    }

    //- receive data sent from processors with greater label
    forAllConstIter(labelHashSet, receiveFrom, iter)
    {
        if( iter.key() <= Pstream::myProcNo() )
            continue;

        List<meshOctreeCubeBasic> mc;

        IPstream fromOtherProc(Pstream::blocking, iter.key());

        fromOtherProc >> mc;

        label currSize = migratedCubes.size();
        migratedCubes.setSize(currSize+mc.size());
        forAll(mc, mcI)
        {
            migratedCubes[currSize] = mc[mcI];
            ++currSize;
        }
    }

    //- send the coordinates of the boxes to processors with lower label
    forAll(sendToProcs, i)
    {
        const label procI = sendToProcs[i];

        if( procI >= Pstream::myProcNo() )
            continue;

        List<meshOctreeCubeBasic> sendCoordinates
        (
            leavesToSend[procI].size()
        );

        forAll(leavesToSend[procI], lI)
        {
            const meshOctreeCube& oc = *leaves[leavesToSend[procI][lI]];
            sendCoordinates[lI] =
                meshOctreeCubeBasic
                (
                    oc.coordinates(),
                    oc.cubeType()
                );
        }

        OPstream toOtherProc
        (
            Pstream::blocking,
            procI,
            sendCoordinates.byteSize()
        );

        toOtherProc << sendCoordinates;
    }

    # ifdef OCTREETiming
    returnReduce(1, sumOp<label>());
    const scalar t4 = omp_get_wtime();
    Info << "Data exchange lasted " << t4-t3 << endl;
    # endif

    //- delete cubes which have been moved to other processors
    octree_.initialCubePtr_->purgeProcessorCubes(Pstream::myProcNo());

    # ifdef OCTREETiming
    returnReduce(1, sumOp<label>());
    const scalar t5 = omp_get_wtime();
    Info << "Purging lasted " << t5-t4 << endl;
    # endif

    //- create boxes from the received coordinates
    forAll(migratedCubes, mcI)
    {
        refineTreeForCoordinates
        (
            migratedCubes[mcI].coordinates(),
            Pstream::myProcNo(),
            migratedCubes[mcI].cubeType()
        );
    }

    createListOfLeaves();

    # ifdef OCTREETiming
    returnReduce(1, sumOp<label>());
    const scalar t6 = omp_get_wtime();
    Info << "Tree refinement lasted " << t6-t5 << endl;
    # endif

    //- update the communication pattern
    updateCommunicationPattern();

    # ifdef OCTREETiming
    returnReduce(1, sumOp<label>());
    const scalar endTime = omp_get_wtime();
    Info << "Updating of communication pattern lasted " << endTime-t6 << endl;
    Info << "Time for load balancing is " << endTime-startTime << endl;
    # endif

    Info << "Finished distributing load between processors" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
