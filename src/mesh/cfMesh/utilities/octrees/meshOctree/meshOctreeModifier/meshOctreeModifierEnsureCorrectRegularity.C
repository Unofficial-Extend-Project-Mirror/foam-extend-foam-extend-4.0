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

void meshOctreeModifier::ensureCorrectRegularity(labelList& refineBox)
{
    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    //- this is needed for parallel runs to reduce the bandwidth
    labelHashSet transferCoordinates;

    labelLongList front;
    forAll(refineBox, leafI)
    {
        if( refineBox[leafI] )
            front.append(leafI);
    }

    DynList<label> neighbours;

    label nMarked;
    do
    {
        nMarked = 0;
        transferCoordinates.clear();
        LongList<meshOctreeCubeCoordinates> processorChecks;

        # ifdef USE_OMP
        # pragma omp parallel private(neighbours)
        # endif
        {
            labelLongList tFront;

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 50)
            # endif
            forAll(front, i)
            {
                const label leafI = front[i];
                const meshOctreeCube* oc = leaves[leafI];

                neighbours.clear();
                octree_.findAllLeafNeighbours(*oc, neighbours);

                forAll(neighbours, neiI)
                {
                    const label nei = neighbours[neiI];

                    if( nei == meshOctreeCubeBasic::OTHERPROC )
                    {
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

                        continue;
                    }

                    if( nei < 0 )
                        continue;

                    const meshOctreeCube& noc = *leaves[nei];

                    if( noc.level() >= oc->level() )
                        continue;

                    if( !refineBox[nei] )
                    {
                        refineBox[nei] = 1;
                        tFront.append(nei);
                    }
                }
            }

            # ifdef USE_OMP

            # pragma omp barrier

            # pragma omp single
            front.clear();

            # pragma omp barrier

            # pragma omp critical
            {
                label start = front.size();
                front.setSize(start+tFront.size());

                forAll(tFront, i)
                    front[start++] = tFront[i];
            }
            # else
            front.transfer(tFront);
            # endif
        }

        nMarked = front.size();

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
            schedule(dynamic, 40) private(neighbours)
            # endif
            forAll(receivedCoords, ccI)
            {
                const meshOctreeCubeCoordinates& cc = receivedCoords[ccI];
                neighbours.clear();
                octree_.findAllLeafNeighbours(cc, neighbours);

                forAll(neighbours, neiI)
                {
                    const label nei = neighbours[neiI];

                    if( nei < 0 )
                        continue;

                    const meshOctreeCube& noc = *leaves[nei];

                    if( noc.level() >= cc.level() )
                        continue;

                    if( !refineBox[nei] )
                    {
                        refineBox[nei] = 1;

                        # ifdef USE_OMP
                        # pragma omp critical
                        # endif
                        front.append(nei);
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

bool meshOctreeModifier::ensureCorrectRegularitySons(labelList& refineBox)
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
