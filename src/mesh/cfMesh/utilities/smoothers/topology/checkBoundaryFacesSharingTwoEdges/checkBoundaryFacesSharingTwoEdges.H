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
    checkBoundaryFacesSharingTwoEdges

Description
    Checks cell connections at the surface. Cells with all vertices at
    the boundary are considered problematic because they get severly skewed
    after being mapped to the geometry surface in most cases. This code performs
    checks of cells with all vertices at the boundary with the rest of the mesh.
    Such cells having a neighbour with all internal faces are removed.

SourceFiles
    checkBoundaryFacesSharingTwoEdges.C

\*---------------------------------------------------------------------------*/

#ifndef checkBoundaryFacesSharingTwoEdges_H
#define checkBoundaryFacesSharingTwoEdges_H

#include "polyMeshGenModifier.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
            Class checkBoundaryFacesSharingTwoEdges Declaration
\*---------------------------------------------------------------------------*/

class checkBoundaryFacesSharingTwoEdges
{
    // Private data
        //- Reference to polyMeshGen
        polyMeshGen& mesh_;

        //- pointer to meshSurfaceEngine
        mutable meshSurfaceEngine* meshSurfacePtr_;

        //- number of boundary faces attached to a boundary vertex
        labelList nBndFacesAtBndPoint_;

        //- a list of boundary points which can be removed from the mesh
        boolList removeBndPoint_;

    // Private member functions
        //- creates meshSurfaceEngine
        void createMeshSurface() const;

        //- returns mesh surface
        inline const meshSurfaceEngine& meshSurface() const
        {
            if( !meshSurfacePtr_ )
                createMeshSurface();

            return *meshSurfacePtr_;
        }

        //- find the faces connected to each boundary edge
        //- mark the points which can be removed
        void findFacesAtBndEdge();

        //- find the number of faces connected to the boundary vertex
        void findBndFacesAtBndVertex();

        //- remove vertices from the faces
        void removeExcessiveVertices();

        //- find boundary faces which shall be decomposed into triangles
        label findBndFacesForDecomposition(boolList& decomposeFace);

        //- Disallow default bitwise copy construct
        checkBoundaryFacesSharingTwoEdges
        (
            const checkBoundaryFacesSharingTwoEdges&
        );

        //- Disallow default bitwise assignment
        void operator=(const checkBoundaryFacesSharingTwoEdges&);

public:

    // Constructors

        //- Construct from polyMeshGen and meshOctree
        checkBoundaryFacesSharingTwoEdges(polyMeshGen& mesh);

    // Destructor

        ~checkBoundaryFacesSharingTwoEdges();

    // Member Functions
        //- find boundary points connected to two boundary faces, only
        void findPoints(labelHashSet& badPoints);

        //- find boundary points connected to two boundary faces
        //- remove points if possible, decompose boundary faces otherwise
        bool improveTopology();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
