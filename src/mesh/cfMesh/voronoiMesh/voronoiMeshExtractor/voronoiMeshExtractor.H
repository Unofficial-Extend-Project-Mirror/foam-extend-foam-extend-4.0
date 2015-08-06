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
    voronoiMeshExtractor

Description
    A class which extracts tet mesh out of an octree structure

SourceFiles
    voronoiMeshExtractor.C

\*---------------------------------------------------------------------------*/

#ifndef voronoiMeshExtractor_H
#define voronoiMeshExtractor_H

#include "polyMeshGenModifier.H"
#include "partTet.H"
#include "tetCreatorOctree.H"
#include "VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class voronoiMeshExtractor Declaration
\*---------------------------------------------------------------------------*/

class voronoiMeshExtractor
{
    // Private data
        //- create tets
        tetCreatorOctree tetCreator_;

        //- reference to the mesh
        polyMeshGen& mesh_;

        //- edges connected to a vertex in the tet mesh
        mutable VRWGraph* pointEdgesPtr_;

        //- edges of the tetrahedrization
        mutable LongList<edge>* edgesPtr_;

        //- tets around each edge
        mutable VRWGraph* edgeTetsPtr_;

        //- is edge at the boundary
        mutable boolList* boundaryEdgePtr_;

    // Private member functions
        //- create and return addressing
        void createAddressing() const;

        const VRWGraph& pointEdges() const;
        const LongList<edge>& edges() const;
        const VRWGraph& edgeTets() const;
        const boolList& boundaryEdge() const;

        //- create points of the voronoi mesh
        void createPoints();

        //- create mesh data
        void createPolyMesh();

        //- delete all pointer data
        void clearOut();

    // Private copy constructor
        //- Disallow default bitwise copy construct
        voronoiMeshExtractor(const voronoiMeshExtractor&);

        //- Disallow default bitwise assignment
        void operator=(const voronoiMeshExtractor&);

    // Private static data
        //- positive orientation
        static label sameOrientation_[6];

        //- opposite orientation
        static label oppositeOrientation_[6];

public:

    // Constructors

        //- Construct from octree and mesh data
        voronoiMeshExtractor
        (
            const meshOctree& octree,
            const IOdictionary& meshDict,
            polyMeshGen& mesh
        );

    // Destructor

        ~voronoiMeshExtractor();


    // Member Functions
        void createMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
