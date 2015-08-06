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

#include "meshOctreeInsideOutside.H"
#include "triSurf.H"
#include "boundBox.H"
#include "labelLongList.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

meshOctreeInsideOutside::meshOctreeInsideOutside
(
    meshOctree& octree
)
:
    octreeModifier_(octree),
    cubeGroup_(octree.numberOfLeaves(), -1),
    cubesInGroup_(),
    groupType_(),
    boundaryDATACubes_(),
    hasOutsideNeighbour_(octree.numberOfLeaves(), false),
    communicationCubes_(),
    neighbouringGroups_()
{
    initialiseBoxes();

    frontalMarking();

    markOutsideCubes();

    reviseDataBoxes();

    markInsideCubes();

    label nInternal(0), nUnknown(0), nData(0), nOutside(0);

    const label nLeaves = octree.numberOfLeaves();
    for(label leafI=0;leafI<nLeaves;++leafI)
    {
        const meshOctreeCubeBasic& oc = octree.returnLeaf(leafI);

        if( oc.cubeType() & meshOctreeCube::INSIDE )
        {
            ++nInternal;
        }
        else if( oc.cubeType() & meshOctreeCube::UNKNOWN )
        {
            ++nUnknown;
        }
        else if( oc.cubeType() & meshOctreeCube::DATA )
        {
            ++nData;
        }
        else if( oc.cubeType() & meshOctreeCube::OUTSIDE )
        {
            ++nOutside;
        }
    }

    if( octree.neiProcs().size() )
    {
        reduce(nInternal, sumOp<label>());
        reduce(nUnknown, sumOp<label>());
        reduce(nData, sumOp<label>());
        reduce(nOutside, sumOp<label>());
    }

    Info << "Number of internal boxes is " << nInternal << endl;
    Info << "Number of outside boxes is " << nOutside << endl;
    Info << "Number of data boxes is " << nData << endl;
    Info << "Number of unknown boxes is " << nUnknown << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

meshOctreeInsideOutside::~meshOctreeInsideOutside()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeInsideOutside::initialiseBoxes()
{
    const LongList<meshOctreeCube*>& leaves = octreeModifier_.leavesAccess();

    # ifdef USE_OMP
    # pragma omp parallel for if( leaves.size() > 1000 )
    # endif
    forAll(leaves, leafI)
    {
        if( leaves[leafI]->hasContainedElements() )
        {
            leaves[leafI]->setCubeType(meshOctreeCubeBasic::DATA);
        }
        else
        {
            leaves[leafI]->setCubeType(meshOctreeCubeBasic::UNKNOWN);
        }
    }
}

void meshOctreeInsideOutside::frontalMarking()
{
    communicationCubes_.clear();
    neighbouringGroups_.clear();

    labelLongList frontCubes;
    DynList<label> neighbours;

    label nGroup(0), nThreads(1);

    const LongList<meshOctreeCube*>& leaves = octreeModifier_.leavesAccess();
    const meshOctree& octree = octreeModifier_.octree();

    boolList commCubes(leaves.size(), false);

    # ifdef USE_OMP
    if( leaves.size() > 1000 )
        nThreads = 3 * omp_get_num_procs();

    # pragma omp parallel num_threads(nThreads) \
    private(frontCubes, neighbours)
    # endif
    {
        LongList<std::pair<label, label> > threadCommPairs;

        # ifdef USE_OMP
        const label threadI = omp_get_thread_num();
        # else
        const label threadI(0);
        # endif

        const label chunkSize = leaves.size() / nThreads + 1;

        const label minLeaf = threadI * chunkSize;

        const label maxLeaf = Foam::min(leaves.size(), minLeaf + chunkSize);

        for(label leafI=minLeaf;leafI<maxLeaf;++leafI)
        {
            if( leaves[leafI]->hasContainedElements() )
                continue;
            if( cubeGroup_[leafI] != -1 )
                continue;

            label groupI;
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            groupI = nGroup++;

            direction cType(meshOctreeCubeBasic::UNKNOWN);
            frontCubes.clear();
            frontCubes.append(leafI);
            cubeGroup_[leafI] = groupI;

            labelLongList neiDATACubes;

            while( frontCubes.size() )
            {
                const label fLabel = frontCubes.removeLastElement();
                octree.findNeighboursForLeaf(fLabel, neighbours);

                forAll(neighbours, neiI)
                {
                    const label nei = neighbours[neiI];
                    if( (nei >= minLeaf) && (nei < maxLeaf) )
                    {
                        if( cubeGroup_[nei] != -1 )
                            continue;

                        if( leaves[nei]->hasContainedElements() )
                        {
                            neiDATACubes.append(nei);
                        }
                        else
                        {
                            frontCubes.append(nei);
                            cubeGroup_[nei] = groupI;
                        }
                    }
                    else if( nei == -1 )
                    {
                        cType = meshOctreeCubeBasic::OUTSIDE;
                    }
                    else if( nei == meshOctreeCubeBasic::OTHERPROC )
                    {
                        commCubes[fLabel] = true;
                    }
                    else
                    {
                        if( leaves[nei]->hasContainedElements() )
                        {
                            neiDATACubes.append(nei);
                        }
                        else
                        {
                            threadCommPairs.append
                            (
                                std::make_pair(fLabel, nei)
                            );
                        }
                    }
                }
            }

            # ifdef USE_OMP
            # pragma omp critical
            # endif
            {
                if( groupI >= boundaryDATACubes_.size() )
                    boundaryDATACubes_.setSize(groupI+1);

                boundaryDATACubes_.setRow(groupI, neiDATACubes);
                groupType_[groupI] = cType;
            }
        }

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        //- find group to neighbouring groups addressing
        List<DynList<label> > localNeiGroups(nGroup);
        forAll(threadCommPairs, cfI)
        {
            const std::pair<label, label>& lp = threadCommPairs[cfI];
            const label groupI = cubeGroup_[lp.first];
            const label neiGroup = cubeGroup_[lp.second];

            if( (neiGroup >= nGroup) || (groupI >= nGroup) )
                FatalError << "neiGroup " << neiGroup
                    << " groupI " << groupI << " are >= than "
                    << "nGroups " << nGroup << abort(FatalError);

            if( neiGroup != -1 )
            {
                localNeiGroups[groupI].appendIfNotIn(neiGroup);
                localNeiGroups[neiGroup].appendIfNotIn(groupI);
            }
        }

        # ifdef USE_OMP
        # pragma omp critical
        # endif
        {
            neighbouringGroups_.setSize(nGroup);

            forAll(localNeiGroups, groupI)
            {
                const DynList<label>& lGroups = localNeiGroups[groupI];

                neighbouringGroups_.appendIfNotIn(groupI, groupI);

                forAll(lGroups, i)
                {
                    neighbouringGroups_.append(groupI, lGroups[i]);
                }
            }
        }
    }

    //- create cubesInGroup_ addressing
    labelList nCubesInGroup(nGroup, 0);
    forAll(cubeGroup_, leafI)
    {
        if( cubeGroup_[leafI] < 0 )
            continue;

        ++nCubesInGroup[cubeGroup_[leafI]];
    }

    cubesInGroup_.setSizeAndRowSize(nCubesInGroup);

    forAllReverse(cubeGroup_, leafI)
    {
        const label groupI = cubeGroup_[leafI];

        if( groupI < 0 )
            continue;

        cubesInGroup_(groupI, --nCubesInGroup[groupI]) = leafI;
    }

    //- mark cubes at inter-processor boundaries
    forAll(commCubes, leafI)
    {
        if( commCubes[leafI] )
            communicationCubes_.append(leafI);
    }

    # ifdef DEBUGSearch
    label nMarked(0);
    forAll(cubeGroup_, leafI)
    {
        if( cubeGroup_[leafI] != -1 )
            ++nMarked;
    }
    reduce(nMarked, sumOp<label>());
    const label totalLeaves = returnReduce(leaves_.size(), sumOp<label>());
    Info << "Total number of leaves " << totalLeaves << endl;
    Info << "Number of marked leaves " << nMarked << endl;
    # endif
}

void meshOctreeInsideOutside::markOutsideCubes()
{
    const LongList<meshOctreeCube*>& leaves = octreeModifier_.leavesAccess();
    const meshOctree& octree = octreeModifier_.octree();

    DynList<label> neighbours;
    label nChanged;
    bool keepUpdating;

    do
    {
        keepUpdating = false;

        do
        {
            nChanged = 0;

            //- make sure that groups created by different threads
            //- have the same information
            forAll(neighbouringGroups_, groupI)
            {
                if( groupType_[groupI] & meshOctreeCubeBasic::OUTSIDE )
                {
                    forAllRow(neighbouringGroups_, groupI, i)
                    {
                        const label neiGroup = neighbouringGroups_(groupI, i);
                        if( groupType_[neiGroup] & meshOctreeCube::UNKNOWN )
                        {
                            ++nChanged;
                            groupType_[neiGroup] = meshOctreeCube::OUTSIDE;
                        }
                    }
                }
            }

            if( nChanged != 0 )
                keepUpdating = true;

        } while( nChanged != 0 );

        do
        {
            nChanged = 0;
            LongList<meshOctreeCubeCoordinates> dataToSend;

            //- go through the list of communicationCubes and send the ones
            //- which are marked as outside
            forAll(communicationCubes_, i)
            {
                const label groupI = cubeGroup_[communicationCubes_[i]];

                if( groupI < 0 )
                    continue;

                if( groupType_[groupI] & meshOctreeCube::OUTSIDE )
                    dataToSend.append(*leaves[communicationCubes_[i]]);
            }

            LongList<meshOctreeCubeCoordinates> receivedCoords;
            octree.exchangeRequestsWithNeighbourProcessors
            (
                dataToSend,
                receivedCoords
            );

            //- go through the list of received coordinates and check if any
            //- local boxes are their neighbours. If a local neighbour is
            //- a DATA box set the hasOutsideNeighbour_ flag to true. If the
            //- local neighbour is of UNKNOWN type set it to OUTSIDE.
            # ifdef USE_OMP
            # pragma omp parallel for if( receivedCoords.size() > 100 ) \
            private(neighbours) schedule(dynamic, 20)
            # endif
            forAll(receivedCoords, i)
            {
                octree.findNeighboursForLeaf(receivedCoords[i], neighbours);

                forAll(neighbours, neiI)
                {
                    const label nei = neighbours[neiI];

                    if( nei < 0 )
                        continue;

                    if( leaves[nei]->hasContainedElements() )
                    {
                        hasOutsideNeighbour_[nei] = true;
                        continue;
                    }
                    if( groupType_[cubeGroup_[nei]] & meshOctreeCube::UNKNOWN )
                    {
                        groupType_[cubeGroup_[nei]] = meshOctreeCube::OUTSIDE;
                        ++nChanged;
                    }
                }
            }

            if( nChanged != 0 )
                keepUpdating = true;

            reduce(nChanged, sumOp<label>());
        } while( nChanged != 0 );

        reduce(keepUpdating, maxOp<bool>());

    } while( keepUpdating );

    //- set OUTSIDE type to the cubes in OUTSIDE groups
    for
    (
        std::map<label, direction>::const_iterator it=groupType_.begin();
        it!=groupType_.end();
        ++it
    )
    {
        if( it->first < 0 )
            continue;

        if( it->second & meshOctreeCubeBasic::OUTSIDE )
        {
            const label groupI = it->first;

            //- set the cube type to OUTSIDE
            forAllRow(cubesInGroup_, groupI, i)
                leaves[cubesInGroup_(groupI, i)]->setCubeType
                (
                    meshOctreeCube::OUTSIDE
                );

            //- set true to the collected DATA boxes
            forAllRow(boundaryDATACubes_, groupI, neiI)
                hasOutsideNeighbour_[boundaryDATACubes_(groupI, neiI)] = true;
        }
    }
}

void meshOctreeInsideOutside::reviseDataBoxes()
{
    //- remove DATA flag from boxes which do not have an OUTSIDE neighbour
    //- and are not surrounded with DATA boxes containing different surface
    //- triangles in different patches
    const LongList<meshOctreeCube*>& leaves = octreeModifier_.leavesAccess();
    const meshOctree& octree = octreeModifier_.octree();
    const triSurf& surface = octree.surface();
    DynList<label> neighbours;

    boolList checkedPatches(leaves.size(), false);

    label nMarked;

    do
    {
        nMarked = 0;

        LongList<meshOctreeCubeCoordinates> checkCoordinates;
        labelHashSet transferCoordinates;

        # ifdef USE_OMP
        # pragma omp parallel for if( leaves.size() > 1000 ) \
        private(neighbours) schedule(dynamic, 20) reduction(+ : nMarked)
        # endif
        forAll(leaves, leafI)
            if( Pstream::parRun() && hasOutsideNeighbour_[leafI] )
            {
                octree.findAllLeafNeighbours(leafI, neighbours);
                forAll(neighbours, neiI)
                    if( neighbours[neiI] == meshOctreeCubeBasic::OTHERPROC )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical
                        # endif
                        {
                            if( !transferCoordinates.found(leafI) )
                            {
                                checkCoordinates.append
                                (
                                    leaves[leafI]->coordinates()
                                );
                                transferCoordinates.insert(leafI);
                            }
                        }

                        break;
                    }
            }
            else if(
                (leaves[leafI]->cubeType() & meshOctreeCube::DATA) &&
                !hasOutsideNeighbour_[leafI]
            )
            {
                meshOctreeCube* oc = leaves[leafI];

                # ifdef DEBUGSearch
                Info << "Box " << leafI << " may not be a DATA box" << endl;
                # endif

                DynList<label> patches;
                const VRWGraph& ct =
                    oc->slotPtr()->containedTriangles_;
                const constRow el = ct[oc->containedElements()];
                forAll(el, elI)
                    patches.appendIfNotIn(surface[el[elI]].region());

                if( patches.size() > 1 )
                    continue;

                checkedPatches[leafI] = true;

                //- check if there exist neighbours
                //- which have some DATA neighbours
                octree.findAllLeafNeighbours(leafI, neighbours);
                forAll(neighbours, neiI)
                {
                    const label nei = neighbours[neiI];

                    if( nei < 0 )
                        continue;

                    if( hasOutsideNeighbour_[nei] )
                    {
                        oc->setCubeType(meshOctreeCube::INSIDE);

                        ++nMarked;
                        break;
                    }
                }
            }

        if( octree.neiProcs().size() )
        {
            LongList<meshOctreeCubeCoordinates> receivedCoords;
            octree.exchangeRequestsWithNeighbourProcessors
            (
                checkCoordinates,
                receivedCoords
            );

            //- check if any of the local neighbours is a data box with
            //- no OUTSIDE neighbours
            # ifdef USE_OMP
            # pragma omp parallel for if( receivedCoords.size() > 100 ) \
            private(neighbours) schedule(dynamic, 20) reduction(+ : nMarked)
            # endif
            forAll(receivedCoords, i)
            {
                octree.findAllLeafNeighbours(receivedCoords[i], neighbours);

                forAll(neighbours, neiI)
                {
                    const label nei = neighbours[neiI];

                    if( nei < 0 )
                        continue;

                    if
                    (
                        (leaves[nei]->cubeType() & meshOctreeCube::DATA) &&
                        !hasOutsideNeighbour_[nei] && checkedPatches[nei]
                    )
                    {
                        leaves[nei]->setCubeType(meshOctreeCube::INSIDE);

                        ++nMarked;
                    }
                }
            }

            reduce(nMarked, sumOp<label>());
        }
    } while( nMarked != 0 );

    # ifdef DEBUGSearch
    label nOutside(0), nData(0), hasOutNei(0);
    forAll(leaves, leafI)
    {
        const direction cType = leaves[leafI]->cubeType();
        if( cType & meshOctreeCubeBasic::OUTSIDE )
            ++nOutside;
        else if( cType & meshOctreeCubeBasic::DATA )
            ++nData;

        if( hasOutsideNeighbour_[leafI] )
            ++hasOutNei;
    }

    reduce(hasOutNei, sumOp<label>());
    reduce(nData, sumOp<label>());
    reduce(nOutside, sumOp<label>());
    Info << "Number of outside boxes " << nOutside << endl;
    Info << "Number of data boxes " << nData << " real data "
        << hasOutNei << endl;
    returnReduce(1, sumOp<label>());
    //::exit(1);
    # endif
}

void meshOctreeInsideOutside::markInsideCubes()
{
    const LongList<meshOctreeCube*>& leaves = octreeModifier_.leavesAccess();
    const meshOctree& octree = octreeModifier_.octree();
    label nChanged;
    bool keepUpdating;
    DynList<label> neighbours;

    //- make INSIDE groups for which it is possible
    for
    (
        std::map<label, direction>::iterator it=groupType_.begin();
        it!=groupType_.end();
        ++it
    )
    {
        const label groupI = it->first;

        if( groupI < 0 )
            continue;

        if( it->second & meshOctreeCubeBasic::UNKNOWN )
        {
            forAllRow(boundaryDATACubes_, groupI, neiI)
            {
                const label cLabel = boundaryDATACubes_(groupI, neiI);
                if(
                    hasOutsideNeighbour_[cLabel] ||
                    (
                        leaves[cLabel]->cubeType() & meshOctreeCube::INSIDE
                    )
                )
                {
                    it->second = meshOctreeCube::INSIDE;
                    break;
                }
            }
        }
    }

    do
    {
        keepUpdating = false;

        //- mark INSIDE groups created by different threads
        do
        {
            nChanged = 0;

            forAll(neighbouringGroups_, groupI)
            {
                if( groupType_[groupI] & meshOctreeCube::INSIDE )
                {
                    forAllRow(neighbouringGroups_, groupI, i)
                    {
                        const label neiGroup = neighbouringGroups_(groupI, i);

                        if( groupType_[neiGroup] & meshOctreeCube::UNKNOWN )
                        {
                            ++nChanged;
                            groupType_[neiGroup] = meshOctreeCube::INSIDE;
                        }
                    }
                }
            }

            if( nChanged != 0 )
                keepUpdating = true;

        } while( nChanged != 0 );

        if( octree.neiProcs().size() == 0 )
            continue;

        //- the code for exchanging data between different processes
        LongList<meshOctreeCubeCoordinates> dataToSend, receivedCoords;

        //- send coordinates of boxes with hasOutsideNeighbour_ flag and
        //- the boxes which have been marked as INSIDE to the neighbouring procs
        forAll(hasOutsideNeighbour_, leafI)
        {
            if(
                hasOutsideNeighbour_[leafI] ||
                (leaves[leafI]->cubeType() & meshOctreeCubeBasic::INSIDE)
            )
            {
                octree.findNeighboursForLeaf(leafI, neighbours);

                forAll(neighbours, neiI)
                    if( neighbours[neiI] == meshOctreeCube::OTHERPROC )
                    {
                        dataToSend.append(leaves[leafI]->coordinates());
                        break;
                    }
            }
        }

        octree.exchangeRequestsWithNeighbourProcessors
        (
            dataToSend,
            receivedCoords
        );

        # ifdef USE_OMP
        # pragma omp parallel for if( receivedCoords.size() > 100 ) \
        private(neighbours) schedule(dynamic, 20)
        # endif
        forAll(receivedCoords, i)
        {
            octree.findNeighboursForLeaf(receivedCoords[i], neighbours);

            forAll(neighbours, neiI)
            {
                const label nei = neighbours[neiI];

                if( nei < 0 )
                    continue;

                const label groupI = cubeGroup_[nei];

                if( groupI < 0 )
                    continue;

                if( groupType_[groupI] & meshOctreeCube::UNKNOWN )
                    groupType_[groupI] = meshOctreeCubeBasic::INSIDE;
            }
        }

        do
        {
            nChanged = 0;
            dataToSend.clear();

            //- go through the list of communicationCubes and send the ones
            //- which are marked as outside
            forAll(communicationCubes_, i)
            {
                if(
                    groupType_[cubeGroup_[communicationCubes_[i]]] &
                    meshOctreeCubeBasic::INSIDE
                )
                    dataToSend.append
                    (
                        leaves[communicationCubes_[i]]->coordinates()
                    );
            }

            receivedCoords.clear();
            octree.exchangeRequestsWithNeighbourProcessors
            (
                dataToSend,
                receivedCoords
            );

            //- go through the list of received coordinates and check if any
            //- local boxes are their neighbours. If a local neighbour is
            //- a DATA box set the hasOutsideNeighbour_ flag to true. If the
            //- local neighbour is of UNKNOWN type set it to OUTSIDE.
            # ifdef USE_OMP
            # pragma omp parallel for if( receivedCoords.size() > 100 ) \
            private(neighbours) schedule(dynamic, 20) reduction(+ : nChanged)
            # endif
            forAll(receivedCoords, i)
            {
                octree.findNeighboursForLeaf(receivedCoords[i], neighbours);

                forAll(neighbours, neiI)
                {
                    const label nei = neighbours[neiI];

                    if( nei < 0 )
                        continue;

                    const label groupI = cubeGroup_[nei];

                    if( groupI < 0 )
                        continue;

                    if( groupType_[groupI] & meshOctreeCube::UNKNOWN )
                    {
                        groupType_[groupI] = meshOctreeCube::INSIDE;
                        ++nChanged;
                    }
                }
            }

            reduce(nChanged, sumOp<label>());

            if( nChanged != 0 )
                keepUpdating = true;

        } while( nChanged != 0 );

        reduce(keepUpdating, maxOp<bool>());

    } while( keepUpdating );

    //- set INSIDE type to the cubes in INSIDE groups
    for
    (
        std::map<label, direction>::const_iterator it=groupType_.begin();
        it!=groupType_.end();
        ++it
    )
    {
        if( it->first < 0 )
            continue;

        if( it->second & meshOctreeCubeBasic::INSIDE )
        {
            const label groupI = it->first;

            //- set the cube type to OUTSIDE
            forAllRow(cubesInGroup_, groupI, i)
                leaves[cubesInGroup_(groupI, i)]->setCubeType
                (
                    meshOctreeCube::INSIDE
                );
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
