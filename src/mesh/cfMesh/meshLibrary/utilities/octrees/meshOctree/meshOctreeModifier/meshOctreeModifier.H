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

Class
    meshOctreeModifier

Description
    Octree for mesh generation

SourceFiles
    meshOctreeModifier.C

\*---------------------------------------------------------------------------*/

#ifndef meshOctreeModifier_H
#define meshOctreeModifier_H

#include "meshOctree.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class meshOctreeModifier Declaration
\*---------------------------------------------------------------------------*/

class meshOctreeModifier
{
    // Private data
        //- Reference to the octree
        meshOctree& octree_;

    // Private member functions
        //- correct refinement such that it produces 1-irregular octree
        void ensureCorrectRegularity(List<direction>& refineBox);

        //- add additional refinement with the following rule. If a son
        //- of a father box is refined, then refine all other sons
        bool ensureCorrectRegularitySons(List<direction>& refineBox);

    // Private copy constructor
        //- Disallow default bitwise copy construct
        meshOctreeModifier(const meshOctreeModifier&);

        //- Disallow default bitwise assignment
        void operator=(const meshOctreeModifier&);

public:

    // Constructors

        //- Construct from octree
        meshOctreeModifier(meshOctree&);

    // Destructor

        ~meshOctreeModifier();


    // Member Functions
        //- return octree
        inline const meshOctree& octree() const;

        //- return neighbour processors
        inline labelList& neiProcsAccess();

        //- return rootBox
        inline boundBox& rootBoxAccess();

        //- return isRootInitialised_
        inline bool& isRootInitialisedAccess();

        //- return search range
        inline scalar& searchRangeAccess();

        //- return the range of leaves at each processor
        inline List<Pair<meshOctreeCubeCoordinates> >& neiRangeAccess();

        //- return initial cube
        inline meshOctreeCube& initialCubeAccess();

        //- return octree slots
        inline List<meshOctreeSlot>& dataSlotsAccess();

        //- return the pointer to the meshOctreeCube at the given position
        inline meshOctreeCube* findCubeForPosition
        (
            const meshOctreeCubeCoordinates&
        ) const;

        //- find leaves contained in a given boundBox
        inline void findLeavesContainedInBox
        (
            const boundBox&,
            DynList<const meshOctreeCube*, 256>&
        ) const;

        //- return leaves
        inline LongList<meshOctreeCube*>& leavesAccess();

        //- create leaves
        inline void createListOfLeaves();

        //- mark additional layers around the leaves selected for refinement
        void markAdditionalLayers
        (
            List<direction>& refineBox,
            const direction nLayers = 1
        ) const;

        //- refine leaves marked for refinement
        //- hexRefinement is activated when it is required to refine all
        //- sons of the same father, if a single son gets marked for refinement
        void refineSelectedBoxes
        (
            List<direction>& refineBox,
            const bool hexRefinement = false
        );

    // functions for parallel runs
        //- distribute leaves of the initial octree to processors
        //- each processor creates a list of neighbouring processors
        void distributeLeavesToProcessors();

        //- move octree cubes from one processor to another
        void loadDistribution(const direction usedType = 0);

        //- refine the tree to add cubes transferred from other processors
        void refineTreeForCoordinates
        (
            const meshOctreeCubeCoordinates&,
            const short procNo = Pstream::myProcNo(),
            const direction cubeType = meshOctreeCubeBasic::UNKNOWN
        );

        void refineTreeForCoordinates
        (
            const meshOctreeCubeCoordinates&,
            const labelList& containedTriangles,
            const labelList& containedEdges,
            const short procNo = Pstream::myProcNo(),
            const direction cubeType = meshOctreeCubeBasic::UNKNOWN

        );

        //- update the communication pattern between the processors
        //- this function is used for load balancing
        void updateCommunicationPattern();

        //- add an additional layer of boxes which belongs to the neighbouring
        //- processors. This simplifies the process of extracting the mesh
        //- template from the octree in a parallel run
        void addLayerFromNeighbouringProcessors();

        //- reduce memory consumption by deleting unnecessary data
        void reduceMemoryConsumption();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshOctreeModifierI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
