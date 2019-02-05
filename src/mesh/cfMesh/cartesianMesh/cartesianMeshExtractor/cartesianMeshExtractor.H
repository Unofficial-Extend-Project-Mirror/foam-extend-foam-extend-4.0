/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

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
