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
    meshSurfaceMapper2D

Description
    Maps vertices of a 2D onto the nearest points on the geometry

SourceFiles
    meshSurfaceMapper2D.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceMapper2D_H
#define meshSurfaceMapper2D_H

#include "labelList.H"
#include "pointField.H"
#include "labelLongList.H"
#include "parMapperHelper.H"
#include "boundBox.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class polyMeshGen;
class polyMeshGen2DEngine;
class meshSurfaceEngine;
class meshSurfacePartitioner;
class triSurfacePartitioner;

/*---------------------------------------------------------------------------*\
                    Class meshSurfaceMapper2D Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceMapper2D
{
    // Private data
        //- reference to mesh surface
        const meshSurfaceEngine& surfaceEngine_;

        //- reference to the octree
        const meshOctree& meshOctree_;

        //- mesh 2D engine
        mutable polyMeshGen2DEngine* mesh2DEnginePtr_;

        //- triSurface partitioner
        mutable triSurfacePartitioner* surfPartitionerPtr_;

        //- mesh surface partitioner
        mutable meshSurfacePartitioner* meshPartitionerPtr_;

        //- labels of edges which are at the boundary which shall be mapped
        //- onto the surface
        labelLongList activeBoundaryEdges_;

    // Private member functions
        //- find active boundary edges
        void findActiveBoundaryEdges();

        //- create a 2D engine
        void create2DEngine() const;

        //- create and return mesh 2D engine
        inline const polyMeshGen2DEngine& mesh2DEngine() const
        {
            if( !mesh2DEnginePtr_ )
                create2DEngine();

            return *mesh2DEnginePtr_;
        }

        //- create surface partitioner
        void createTriSurfacePartitioner() const;

        //- create and reurn surface partitioner
        inline const triSurfacePartitioner& surfacePartitioner() const
        {
            if( !surfPartitionerPtr_ )
                createTriSurfacePartitioner();

            return *surfPartitionerPtr_;
        }

        //- create mesh surface partitioner
        void createMeshSurfacePartitioner() const;

        //- create and return mesh surface partitioner
        inline const meshSurfacePartitioner& meshPartitioner() const
        {
            if( !meshPartitionerPtr_ )
                createMeshSurfacePartitioner();

            return *meshPartitionerPtr_;
        }

        //- delete dynamically allocated data
        void clearOut();

    // Private member functions needed for parallel execution
        //- find mapping distance for corner points
        void findMappingDistance
        (
            const labelLongList&,
            std::map<label, scalar>&
        ) const;

        //- map to the location with the smallest distance
        void mapToSmallestDistance(LongList<parMapperHelper>&);

        //- Disallow default bitwise copy construct
        meshSurfaceMapper2D(const meshSurfaceMapper2D&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceMapper2D&);

public:

    // Constructors

        //- Construct from meshSurfaceEngine and octree
        meshSurfaceMapper2D(const meshSurfaceEngine&, const meshOctree&);

    // Destructor

        ~meshSurfaceMapper2D();

    // Member Functions

        //- adjust z coordinates of the mesh to the ones in the surface mesh
        void adjustZCoordinates();

        //- projects surface vertices onto their nearest location
        //- on the surface mesh
        void mapVerticesOntoSurface();

        //- projects selected edge to their nearest location
        //- on the surface mesh
        void mapVerticesOntoSurface(const labelLongList& edgesToMap);

        //- projects edges with corners onto their nearest location
        //- on the surface mesh
        void mapCorners();

        //- projects selected edges with corners onto their nearest
        //- locations on the surface mesh
        void mapCorners(const labelLongList& edgesToMap);

        //- projects surface edges onto the surface with respect
        //- to the surface patch they belong to. Corners are respected
        void mapVerticesOntoSurfacePatches();

        //- projects selected surface edges onto the surface with respect
        //- to the surface patch they belong to. Corners are respected
        void mapVerticesOntoSurfacePatches(const labelLongList& edgesToMap);

        //- a combination of mapping and smoothing intended for better
        //- feature capturing
        void preMapVertices(const label nIterations = 2);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
