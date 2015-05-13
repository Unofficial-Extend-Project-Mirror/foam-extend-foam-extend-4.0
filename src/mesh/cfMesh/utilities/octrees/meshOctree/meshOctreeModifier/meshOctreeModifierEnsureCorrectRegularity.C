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
#include "HashSet.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeModifier::ensureCorrectRegularity(List<direction>& refineBox)
{
    const FixedList<meshOctreeCubeCoordinates, 26>& rp =
        octree_.regularityPositions_;
    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    //- this is needed for parallel runs to reduce the bandwidth
    labelHashSet transferCoordinates;

    labelLongList front;
    forAll(refineBox, leafI)
    {
        if( refineBox[leafI] )
            front.append(leafI);
    }

    FixedList<meshOctreeCube*, 26> neighbours;

    label nMarked;
    do
    {
        nMarked = 0;
        LongList<meshOctreeCubeCoordinates> processorChecks;

        # ifdef USE_OMP
        # pragma omp parallel if( front.size() > 1000 ) \
        private(neighbours) reduction(+ : nMarked)
        # endif
        {
            labelLongList tFront;

            # ifdef USE_OMP
            # pragma omp for
            # endif
            forAll(front, i)
                tFront.append(front[i]);

            # ifdef USE_OMP
            # pragma omp barrier
            # endif

            front.clear();

            while( tFront.size() != 0 )
            {
                const label leafI = tFront.removeLastElement();
                const meshOctreeCube* oc = leaves[leafI];

                forAll(rp, posI)
                {
                    const meshOctreeCubeCoordinates cc
                    (
                        oc->coordinates() + rp[posI]
                    );

                    const label neiLabel = octree_.findLeafLabelForPosition(cc);

                    if( neiLabel > -1 )
                    {
                        neighbours[posI] = leaves[neiLabel];
                    }
                    else if( neiLabel == -1 )
                    {
                        neighbours[posI] = NULL;
                    }
                    else if( neiLabel == meshOctreeCubeBasic::OTHERPROC )
                    {
                        neighbours[posI] = NULL;

                        # ifdef USE_OMP
                        # pragma omp critical
                        # endif
                        {
                            if( !transferCoordinates.found(leafI) )
                            {
                                processorChecks.append(oc->coordinates());
                                transferCoordinates.insert(leafI);
                            }
                        }
                    }
                }

                forAll(neighbours, neiI)
                {
                    const meshOctreeCube* nei = neighbours[neiI];
                    if( !nei ) continue;
                    if( !nei->isLeaf() ) continue;
                    if( nei->level() >= oc->level() ) continue;

                    if( !refineBox[nei->cubeLabel()] )
                    {
                        refineBox[nei->cubeLabel()] = 1;
                        tFront.append(nei->cubeLabel());
                    }
                }
            }
        }

        if( octree_.neiProcs().size() )
        {
            LongList<meshOctreeCubeCoordinates> receivedCoords;
            octree_.exchangeRequestsWithNeighbourProcessors
            (
                processorChecks,
                receivedCoords
            );

            //- check consistency with received cube coordinates
            # ifdef USE_OMP
            # pragma omp parallel for if( receivedCoords.size() > 100 ) \
            schedule(dynamic, 40)
            # endif
            forAll(receivedCoords, ccI)
            {
                forAll(rp, posI)
                {
                    const meshOctreeCubeCoordinates cc
                    (
                        receivedCoords[ccI] + rp[posI]
                    );

                    const meshOctreeCube* nei =
                        octree_.findCubeForPosition(cc);

                    if( !nei ) continue;
                    if( !nei->isLeaf() ) continue;
                    if( nei->level() >= cc.level() ) continue;

                    if( !refineBox[nei->cubeLabel()] )
                    {
                        refineBox[nei->cubeLabel()] = 1;

                        # ifdef USE_OMP
                        # pragma omp critical
                        # endif
                        front.append(nei->cubeLabel());
                    }
                }
            }

            nMarked = front.size();

            //- calculate the number of selected boxes over all processors
            reduce(nMarked, sumOp<label>());
        }
    }
    while( nMarked != 0 );
}

bool meshOctreeModifier::ensureCorrectRegularitySons(List<direction>& refineBox)
{
    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    LongList<meshOctreeCubeCoordinates> transferCoordinates;

    label nMarked(0);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100) reduction(+ : nMarked)
    # endif
    forAll(leaves, leafI)
    {
        if( !refineBox[leafI] )
            continue;

        const meshOctreeCubeCoordinates cc = leaves[leafI]->reduceLevelBy(1);

        for(label scI=0;scI<8;++scI)
        {
            const label neiLeaf =
                octree_.findLeafLabelForPosition(cc.refineForPosition(scI));

            if( neiLeaf >= 0 && !refineBox[neiLeaf] )
            {
                //- mark this leaf for refinement
                ++nMarked;
                refineBox[neiLeaf] = 1;
            }
            else if( neiLeaf == meshOctreeCube::OTHERPROC )
            {
                //- propagate this information to other processors
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                transferCoordinates.append(cc);
            }
        }
    }

    if( octree_.neiProcs().size() )
    {
        LongList<meshOctreeCubeCoordinates> receivedCoords;
        octree_.exchangeRequestsWithNeighbourProcessors
        (
            transferCoordinates,
            receivedCoords
        );

        # ifdef USE_OMP
        # pragma omp parallel for if( receivedCoords.size() > 100 ) \
        reduction(+ : nMarked)
        # endif
        forAll(receivedCoords, ccI)
        {
            const meshOctreeCubeCoordinates& cc = receivedCoords[ccI];

            for(label scI=0;scI<8;++scI)
            {
                const label neiLeaf =
                    octree_.findLeafLabelForPosition(cc.refineForPosition(scI));

                if( neiLeaf >= 0 && !refineBox[neiLeaf] )
                {
                    //- mark this leaf for refinement
                    ++nMarked;
                    refineBox[neiLeaf] = 1;
                }
            }
        }
    }

    reduce(nMarked, sumOp<label>());

    if( nMarked )
        return true;

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
