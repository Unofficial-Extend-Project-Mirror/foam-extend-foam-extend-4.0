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
#include "triSurf.H"
#include "HashSet.H"

# ifdef USE_OMP
#include <omp.h>
# endif

#include <sys/stat.h>

//#define OCTREETiming
//#define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeModifier::markAdditionalLayers
(
    List<direction>& refineBox,
    const direction nLayers
) const
{
    const FixedList<meshOctreeCubeCoordinates, 26>& rp =
        octree_.regularityPositions_;
    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    //- this is needed for parallel runs to reduce the bandwidth
    labelHashSet transferCoordinates;

    FixedList<meshOctreeCube*, 26> neighbours;

    for(direction i=1;i<=nLayers;++i)
    {
        LongList<meshOctreeCubeCoordinates> processorChecks;

        # ifdef USE_OMP
        # pragma omp parallel for if( leaves.size() > 1000 ) \
        private(neighbours) schedule(dynamic, 20)
        # endif
        forAll(leaves, leafI)
        {
            if( refineBox[leafI] != i )
                continue;

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
                if( nei->level() > oc->level() ) continue;

                if( !refineBox[nei->cubeLabel()] )
                {
                    refineBox[nei->cubeLabel()] = i+1;
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
            # pragma omp parallel for if( receivedCoords.size() > 1000 ) \
            schedule(dynamic, 20)
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
                    if( nei->level() > cc.level() ) continue;

                    if( !refineBox[nei->cubeLabel()] )
                    {
                        refineBox[nei->cubeLabel()] = i+1;
                    }
                }
            }
        }
    }
}

void meshOctreeModifier::refineSelectedBoxes
(
    List<direction>& refineBox,
    const bool hexRefinement
)
{
    # ifdef OCTREETiming
    const scalar startTime = omp_get_wtime();
    # endif

    //- ensure that refinement will produce 1-irregular octree
    do
    {
        ensureCorrectRegularity(refineBox);
    } while( hexRefinement && ensureCorrectRegularitySons(refineBox) );

    # ifdef OCTREETiming
    const scalar regTime = omp_get_wtime();
    Info << "Time for ensuring regularity " << (regTime-startTime) << endl;
    # endif

    const triSurf& surface = octree_.surface();
    const boundBox& rootBox = octree_.rootBox();
    const LongList<meshOctreeCube*>& leaves = octree_.leaves_;

    //- this is needed for thread safety
    //- such solutions make me a sad bunny :(
    surface.facetEdges();
    surface.edgeFacets();
    surface.edges();

    # ifdef USE_OMP
    # pragma omp parallel num_threads(octree_.dataSlots_.size())
    # endif
    {
        # ifdef USE_OMP
        meshOctreeSlot* slotPtr = &octree_.dataSlots_[omp_get_thread_num()];
        # else
        meshOctreeSlot* slotPtr = &octree_.dataSlots_[0];
        # endif

        if( !octree_.isQuadtree() )
        {
            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 100)
            # endif
            forAll(leaves, leafI)
            {
                if( refineBox[leafI] )
                    leaves[leafI]->refineCube(surface, rootBox, slotPtr);
            }
        }
        else
        {
            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 100)
            # endif
            forAll(leaves, leafI)
            {
                if( refineBox[leafI] )
                    leaves[leafI]->refineCube2D(surface, rootBox, slotPtr);
            }
        }
    }

    createListOfLeaves();

    # ifdef OCTREETiming
    Info << "Time for actual refinement " << (omp_get_wtime()-regTime) << endl;
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
