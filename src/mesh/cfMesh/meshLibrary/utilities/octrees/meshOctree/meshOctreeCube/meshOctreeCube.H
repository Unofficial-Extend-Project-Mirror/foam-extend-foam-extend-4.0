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
    meshOctreeCube

Description
    A cube stores information needed for mesh generation

SourceFiles
    meshOctreeCube.C

\*---------------------------------------------------------------------------*/

#ifndef meshOctreeCube_H
#define meshOctreeCube_H

#include "DynList.H"
#include "LongList.H"
#include "boolList.H"
#include "meshOctreeCubeBasic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Ostream;
class triSurf;
class VRWGraph;
class meshOctreeSlot;

/*---------------------------------------------------------------------------*\
                    Class meshOctreeCube Declaration
\*---------------------------------------------------------------------------*/

class meshOctreeCube
: public meshOctreeCubeBasic
{
    // Private data
        //- pointer the slot containing this cube
        meshOctreeSlot* activeSlotPtr_;

        //- pointer to the first child element
        meshOctreeCube** subCubesPtr_;

        //- position of the cube in the list of leaves
        mutable label cubeLabel_;

        //- label of the row which contains elements contained in the cube
        label containedElementsLabel_;

        //- labels of contained surface edges
        label containedEdgesLabel_;

    // Static data
        //- this data is needed to map the the current refinement pattern
        //- to the Hilbert space-filling curve in order to reduce the bandwidth
        //- in the list of leaves. The procedure is performed using the table
        //- published in (Campbell etal: Dynamic Octree Load Balancing Using
        //- Space-Filling Curves) which is adjusted to the current ref pattern
        static const label hOrder_[24][8];

        //- orientation of the Hilbert curve in the refined cubes of a cube
        //- with the ordering in the above table
        static const label hOrient_[24][8];

    // Private member functions
        //- find edges contained in the cube
        void findContainedEdges
        (
            const triSurf&,
            const boundBox&
        );

        //- store data which is contained in the child cube
        //- this is needed for parallel coarsening of the octree
        //inline void reclaimDataFromChild(const label scI);

        //- Disallow default bitwise copy construct
        meshOctreeCube(const meshOctreeCube&);

        //- Disallow copy construct from meshOctreeCubeBasic
        meshOctreeCube(const meshOctreeCubeBasic&);

public:

    // Constructors

        //- Default constructor
        inline meshOctreeCube();

        //- Construct from coordinates
        meshOctreeCube(const meshOctreeCubeCoordinates&);

        //- Construct from coordinates and the number of surface triangles
        //- this constructor is used for the initial cube, only
        meshOctreeCube
        (
            const meshOctreeCubeCoordinates& cc,
            const label nElmts,
            meshOctreeSlot* slotPtr
        );

    // Destructor

        ~meshOctreeCube();

    // Member Functions

        //- refine cube in two directions, it is used for generating quadtrees
        void refineCube2D
        (
            const triSurf&,
            const boundBox&,
            meshOctreeSlot* slotPtr = NULL
        );

        //- subdivide the octree cube
        void refineCube
        (
            const triSurf&,
            const boundBox&,
            meshOctreeSlot* slotPtr = NULL
        );

        //- create missing child cubes and refine them until the required
        //- cube is reached (this function is used for parallel octree creation)
        void refineMissingCube
        (
            const triSurf&,
            const boundBox&,
            const label scI,
            meshOctreeSlot* slotPtr = NULL
        );

        //- create missing child cubes and refine them until the required
        //- cube is reached (this function is used for parallel octree creation)
        //- contained elements and edges are provided as argument
        void refineMissingCube
        (
            const label scI,
            const label elementsRowI = -1,
            const label edgesRowI = -1,
            meshOctreeSlot* slotPtr = NULL
        );

        //- return the pointer to the slot containing the cube
        inline const meshOctreeSlot* slotPtr() const;

        //- position of the cube in the list of leaves
        inline label cubeLabel() const;

        //- check if the cube is a leaf
        inline bool isLeaf() const;

        //- return a pointer to a child cubes at given position
        inline meshOctreeCube* subCube(const label) const;

        //- return the pointers to the sons
        FixedList<meshOctreeCube*, 8> subCubes() const;

        //- check if this box has some contained triangles
        bool hasContainedTriangles
        (
            const triSurf&,
            const boundBox&,
            const VRWGraph& containedElements
        ) const;

        //- return true if the box contains some triangles
        inline bool hasContainedElements() const;

        //- returns the list of contained surface elements
        //- this list is only allocated for cubes containing elements
        inline label containedElements() const;

        //- returns true if the box contains some edges
        inline bool hasContainedEdges() const;

        //- returns the list of contained surface edges
        //- this list is only allocated for cubes containing edges
        inline label containedEdges() const;

    // Recursive functions
        //- leaves contained in the given box
        void leavesInBox
        (
            const boundBox& rootBox,
            const boundBox& searchingBox,
            DynList<const meshOctreeCube*, 256>&
        ) const;

        //- find leaves for a given cube
        void findLeaves
        (
            LongList<meshOctreeCube*>& leaves
        ) const;

        //- find coordinates of cubes which are located on other processors
        void findCoordinatesOfMissingCubes
        (
            LongList<meshOctreeCubeCoordinates>& coordinates
        ) const;

        //- count number of originating from this cube
        void countChildCubes(label& nCubes) const;

        //- delete boxes which are not local to the given processor
        bool purgeProcessorCubes(const short procNo);

    // Member operators

        //- assignment
        inline void operator=(const meshOctreeCube&);
        friend Ostream& operator<<(Ostream&, const meshOctreeCube&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshOctreeCubeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
