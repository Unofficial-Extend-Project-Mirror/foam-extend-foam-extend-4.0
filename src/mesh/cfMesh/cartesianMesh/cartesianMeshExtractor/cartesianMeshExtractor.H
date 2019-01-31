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
    cartesianMeshExtractor

Description
    Generates the cartesian mesh template from the octree

SourceFiles
    cartesianMeshExtractor.C

\*---------------------------------------------------------------------------*/

#ifndef cartesianMeshExtractor_H
#define cartesianMeshExtractor_H

#include "polyMeshGenModifier.H"
#include "meshOctreeAddressing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
class IOdictionary;

/*---------------------------------------------------------------------------*\
                    Class cartesianMeshExtractor Declaration
\*---------------------------------------------------------------------------*/

class cartesianMeshExtractor
{
    // Private data
        //- reference to the octree addressing
        meshOctreeAddressing octreeCheck_;
    
        //- reference to the mesh
        polyMeshGen& mesh_;
    
        //- decompose split hex cells
        bool decomposeSplitHexes_;
    
        //- cell label for a given leaf
        labelList* leafCellLabelPtr_;

    // Private member functions
        //- delete all freestore data
        void clearOut();

        //- create vertices and pointLeaves addressing
        void createPointsAndAddressing();

        //- create mesh data
        void createPolyMesh();
        
        //- decompose split hexes into pyramids and tets
        void decomposeSplitHexesIntoTetsAndPyramids();

    // Private copy constructor
        //- Disallow default bitwise copy construct
        cartesianMeshExtractor(const cartesianMeshExtractor&);

        //- Disallow default bitwise assignment
        void operator=(const cartesianMeshExtractor&);

public:

    // Constructors

        //- Construct from octree and mesh data
        cartesianMeshExtractor
        (
            meshOctree& octree,
            const IOdictionary& meshDict,
            polyMeshGen& mesh
        );

    // Destructor

        ~cartesianMeshExtractor();


    // Member Functions

        //- decompose split hexes into standard cells
        void decomposeSplitHexes();
        
        //- create the mesh with the above options
        void createMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
