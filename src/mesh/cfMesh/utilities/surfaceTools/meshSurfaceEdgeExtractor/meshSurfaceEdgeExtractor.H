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
    meshSurfaceEdgeExtractor

Description
    Extracts edges and corners to cornform the mesh to the geometry

SourceFiles
    meshSurfaceEdgeExtractor.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceEdgeExtractor_H
#define meshSurfaceEdgeExtractor_H

#include "polyMeshGenModifier.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class meshSurfaceEdgeExtractor Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceEdgeExtractor
{
    // Private data
        //- mesh
        polyMeshGen& mesh_;
    
        label nPoints_;
    
        boolList boundaryCell_;
    
        List<direction> nFacesInCell_;

        //- reference to the octree
        const meshOctree& meshOctree_;

        //- regions for boundary vertices
        VRWGraph pointRegions_;

    // Private member functions
        //- remove old boundary faces
        void removeOldBoundaryFaces();
    
        //- create vertices on surface edges
        void createEdgeVertices();
    
        //- create new boundary faces
        void createBoundaryFaces();

        //- Disallow default bitwise copy construct
        meshSurfaceEdgeExtractor(const meshSurfaceEdgeExtractor&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceEdgeExtractor&);

public:

    // Constructors

        //- Construct from mesh data
        meshSurfaceEdgeExtractor
        (
            polyMeshGen& mesh,
            const meshOctree& octree,
            const labelList& pointRegion
        );

    // Destructor

        ~meshSurfaceEdgeExtractor();

    // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
