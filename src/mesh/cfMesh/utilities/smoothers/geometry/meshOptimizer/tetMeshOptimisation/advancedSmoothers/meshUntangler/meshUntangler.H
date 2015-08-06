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
    meshUntangler

Description
    Mesh smoothing without any topological changes. The vertex is put
    into the centre of the convex feasible region.

SourceFiles
    meshUntangler.C

\*---------------------------------------------------------------------------*/

#ifndef meshUntangler_H
#define meshUntangler_H

#include "simplexSmoother.H"
#include "polyMeshGen.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class plane;
class partTetMeshSimplex;

/*---------------------------------------------------------------------------*\
                        Class meshUntangler Declaration
\*---------------------------------------------------------------------------*/
        
class meshUntangler
:
    public simplexSmoother
{
    // Private nested classes
    
        class cutRegion
        {
            // Private data
                DynList<point, 64>* pointsPtr_;
                DynList<edge, 128>* edgesPtr_;
                DynList<DynList<label, 8>, 64>* facesPtr_;
            
                //- helper data
                DynList<point, 64>* cPtsPtr_;
                DynList<edge, 128>* cEdgesPtr_;
                DynList<DynList<label, 8>, 64>* cFacesPtr_;
            
                DynList<label, 64> newVertexLabel_;
                DynList<scalar, 64> vertexDistance_;
                DynList<direction, 64> vertexTypes_;
                DynList<label, 128> newEdgeLabel_;
                label origNumVertices_;
            
                scalar tol_;
                bool valid_;
            
            // Private member functions
                bool findNewVertices(const plane& plane);
                void findNewEdges();
                void findNewFaces();
            
                //- tie break to resolve problems due to truncation
                //- errors
                void tieBreak(const DynList<label, 8>& f);
                
                //- remove coincident vertices to improve tie breaking
                void removeCoincidentVertices();
            
                void createInitialConfiguration(const boundBox&);
                
                // Enumerations
                enum vTypes
                {
                    NONE = 0,
                    KEEP = 1,
                    INPLANE = 2,
                    USED = 4
                };
            
            public:
                
            // Constructor
                //- Construct from boundBox
                cutRegion(const boundBox&);
            
            // Destructor
                ~cutRegion();
            
            // Member functions
                //- return the vertices of the feasible region
                inline const DynList<point, 64>& points() const
                {
                    return *pointsPtr_;
                }
                
                //- cut the region woth the plane
                void planeCut(const plane& plane);
                
                //- export the feasible region as polyMeshGen
                void createPolyMeshFromRegion(polyMeshGen&) const;
        };
    
    public:
        
    // Constructor
        //- Construct from partTetMeshSimplex
        meshUntangler(partTetMeshSimplex& simplex);
    
    // Destructor
        ~meshUntangler();
    
    // Member functions
        //- improve the position the the node by putting it
        //- into the feasible region
        void optimizeNodePosition(const scalar tol = 0.001);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
