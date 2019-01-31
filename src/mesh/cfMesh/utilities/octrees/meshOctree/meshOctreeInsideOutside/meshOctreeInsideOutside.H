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
    meshOctreeInsideOutside
Description
    Creates inside/outside information for the octree

SourceFiles
    meshOctreeInsideOutside.C

\*---------------------------------------------------------------------------*/

#ifndef meshOctreeInsideOutside_H
#define meshOctreeInsideOutside_H

#include "meshOctreeModifier.H"
#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class meshOctreeInsideOutside Declaration
\*---------------------------------------------------------------------------*/

class meshOctreeInsideOutside
{
    // Private data
    
        //- meshOctreeAddressing
        meshOctreeModifier octreeModifier_;
    
        //- group for a given leaf
        labelLongList cubeGroup_;
    
        //- cubes belonging to each group of octree boxes
        VRWGraph cubesInGroup_;
    
        //- type for a group
        std::map<label, direction> groupType_;
    
        //- boundary DATA boxes for a given group
        VRWGraph boundaryDATACubes_;
    
        //- flag for DATA boxes next to the OUTSIDE boxes
        boolList hasOutsideNeighbour_;
    
        //- label of cubes at processor boundaries
        labelLongList communicationCubes_;
    
        //- labels of cubes marked by different threads
        VRWGraph neighbouringGroups_;
    
    // Private member functions
    
        //- initialise octree boxes
        void initialiseBoxes();
    
        //- frontal marking algorithm
        void frontalMarking();
        
        //- mark OUTSIDE octree boxes
        void markOutsideCubes();
        
        //- revise DATA boxes
        void reviseDataBoxes();
        
        //- mark INSIDE octree boxes
        void markInsideCubes();
    
    public:
        
    // Constructor
    
        //- Construct from octree
        meshOctreeInsideOutside
        (
            meshOctree& octree
        );
    
    // Destructor
    
        ~meshOctreeInsideOutside();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
