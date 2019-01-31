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
    checkIrregularSurfaceConnections

Description
    Check mesh surface for for connections which cannot exist in a valid mesh.
    In case some irregular connections exist, it is possible to fix them.

SourceFiles
    checkIrregularSurfaceConnections.C

\*---------------------------------------------------------------------------*/

#ifndef checkIrregularSurfaceConnections_H
#define checkIrregularSurfaceConnections_H

#include "polyMeshGenModifier.H"
#include "meshSurfaceEngine.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class checkIrregularSurfaceConnections Declaration
\*---------------------------------------------------------------------------*/

class checkIrregularSurfaceConnections
{
    // Private data
        //- Reference to polyMeshGen
        polyMeshGen& mesh_;
    
        //- pointer to meshSurfaceEngine
        mutable meshSurfaceEngine* meshSurfacePtr_;

    // Private member functions
        //- create meshSurfaceEngine
        inline const meshSurfaceEngine& surfaceEngine() const
        {
            if( !meshSurfacePtr_ )
                meshSurfacePtr_ = new meshSurfaceEngine(mesh_);
            
            return *meshSurfacePtr_;
        }
        
        //- clear meshSurfaceEngine. Used in case if the mesh is modified
        inline void clearMeshEngine()
        {
            deleteDemandDrivenData(meshSurfacePtr_);
        }
        
        //- check if there exist vertices with more than one groups of cells
        //- attached to it
        bool checkAndFixCellGroupsAtBndVertices
        (
            labelHashSet& badVertices,
            const bool removeConnections = false
        );
        
        //- check if all boundary edges are connected
        //- to two boundary faces, only!
        bool checkEdgeFaceConnections
        (
            labelHashSet& badVertices,
            const bool removeCells = false
        );
        
        //- check if all boundary face connected to a boundary vertex
        //- are connected into a single group of faces that can be visited
        //- over edges
        bool checkFaceGroupsAtBndVertices
        (
            labelHashSet& badVertices,
            const bool removeCells = false
        );

        //- Disallow default bitwise copy construct
        checkIrregularSurfaceConnections
        (
            const checkIrregularSurfaceConnections&
        );

        //- Disallow default bitwise assignment
        void operator=(const checkIrregularSurfaceConnections&);

public:

    // Constructors

        //- Construct from polyMeshGen
        checkIrregularSurfaceConnections(polyMeshGen& mesh);

    // Destructor

        ~checkIrregularSurfaceConnections();

    // Member Functions
        //- find vertices where the surface is not valid
        void checkIrregularVertices(labelHashSet& badVertices);

        //- find invalid connections and fix them on the fly
        bool checkAndFixIrregularConnections();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
