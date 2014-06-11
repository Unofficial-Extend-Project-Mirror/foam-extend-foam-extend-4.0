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
    meshOptimizer

Description
    Mesh smoothing without any topological changes

SourceFiles
    meshOptimizer.C

\*---------------------------------------------------------------------------*/

#ifndef meshOptimizer_H
#define meshOptimizer_H

#include "polyMeshGen.H"
#include "boolList.H"
#include "VRWGraph.H"
#include "DynList.H"
#include "partTet.H"
#include "HashSet.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;
class plane;
class partTetMesh;

/*---------------------------------------------------------------------------*\
                        Class meshOptimizer Declaration
\*---------------------------------------------------------------------------*/

class meshOptimizer
{
    // Private data
        //- reference to the mesh
        polyMeshGen& mesh_;
    
        //- location of vertex (internal, boundary, edge, corner)
        List<direction> vertexLocation_;
    
        //- mesh surface
        mutable meshSurfaceEngine* msePtr_;

    // Private member functions
        //- return mesh surface
        const meshSurfaceEngine& meshSurface() const;
        void clearSurface();
    
        //- find problematic faces
        label findBadFaces(labelHashSet&, const boolList&) const;
        label findLowQualityFaces(labelHashSet&, const boolList&) const;
    
    // Nested classes
    
        class laplaceSmoother
        {
            // Private data
                //- reference to the mesh
                polyMeshGen& mesh_;
            
                //- location of vertex (internal, boundary, edge, corner)
                const List<direction>& vertexLocation_;
            
            // Private member functions
                //- smooth the node using the laplacian smoother
                //- new position is the average of the neighbouring vertices
                void laplacian(const labelLongList&, const label);
                void laplacianSurface(const labelLongList&, const label);
            
                void laplacianParallel
                (
                    const labelLongList& procPoints,
                    const bool smoothOnlySurfaceNodes = false
                );
                
                //- smooth the node using the laplacian smoother
                //- new position is the average of the centres of faces attached
                //- to the vertex
                void laplacianPC(const labelLongList&, const label);
                void laplacianPCParallel(const labelLongList& procPoints);
                
                //- smooth the node using the laplacian smoother
                //- new position is the average of the centres of faces attached
                //- to the vertex
                void laplacianWPC(const labelLongList&, const label);
                void laplacianWPCParallel(const labelLongList& procPoints);
                
                //- update geometry after smoothing
                void updateMeshGeometry(const labelLongList& smoothPoints);
                
                //- Disallow default bitwise copy construct
                laplaceSmoother(const laplaceSmoother&);

                //- Disallow default bitwise assignment
                void operator=(const laplaceSmoother&);
    
        public:
                
            // Constructor
        
                //- Construct from mesh and vertex locations
                laplaceSmoother(polyMeshGen&, const List<direction>&);
        
            // Destructor
        
                ~laplaceSmoother();
        
            // Member Functions
                //- new position is the average of the neighbouring vertices
                void optimizeLaplacian(const label nIterations = 1);
                void optimizeLaplacian
                (
                    const labelHashSet& badFaces,
                    const label nIterations = 1
                );
        
                //- new position of surface point is the average of
                //- the neighbouring surface vertices
                void optimizeSurfaceLaplacian
                (
                    const labelHashSet& badFaces,
                    const label nIterations = 1
                );
                
                //- new positions are the average of the centres of the cells
                //- adjacent to the vertex
                void optimizeLaplacianPC(const label nIterations = 1);
                void optimizeLaplacianPC
                (
                    const labelHashSet& badFaces,
                    const label nIterations = 1
                );
                
                //- new positions are the average of the centres of the cells
                //- adjacent to the vertex weighted by cell volumes
                void optimizeLaplacianWPC(const label nIterations = 1);
                void optimizeLaplacianWPC
                (
                    const labelHashSet& badFaces,
                    const label nIterations = 1
                );
        };
        
        //- Disallow default bitwise copy construct
        meshOptimizer(const meshOptimizer&);

        //- Disallow default bitwise assignment
        void operator=(const meshOptimizer&);
        
    // enumerators
        enum vertexType_
        {
            INSIDE = 1,
            BOUNDARY = 2,
            EDGE = 4,
            CORNER = 8,
            PARALLELBOUNDARY = 16
        };

public:

    // Constructors

        //- Construct from mesh
        meshOptimizer(polyMeshGen& mesh);


    // Destructor

        ~meshOptimizer();

    // Member Functions
        //- smooth surface vertices
        void optimizeSurface(const meshOctree&);

        //- performs mesh untangling based on detected negative normals
        void untangleMeshFV();

        //- performs mesh optimisation for faces with high non-orthogonality
        //- and skewness
        void optimizeLowQualityFaces();

        //- final optimisation
        void optimizeMeshFV();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
