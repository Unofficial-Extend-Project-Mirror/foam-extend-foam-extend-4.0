/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::walkPatch

Description
    Collection of static functions to do various simple patch related things.

SourceFiles
    walkPatch.C

\*---------------------------------------------------------------------------*/

#ifndef walkPatch_H
#define walkPatch_H

#include "labelList.H"
#include "primitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class walkPatch Declaration
\*---------------------------------------------------------------------------*/

class walkPatch
{
    // Private data

        //- Reference to patch to walk on
        const primitivePatch& pp_;

        //- Reference to zones
        const labelList& faceZone_;

        //- How to walk through faces
        const bool reverse_;

        //- Reference to list to mark off visited faces
        boolList& visited_;


        // Faces visited
        dynamicLabelList visitOrder_;

        // Index in face of vertex it was visited through
        dynamicLabelList indexInFace_;


    // Private Member Functions

        //- Get other face using v0, v1. Returns -1 if none.
        label getNeighbour
        (
            const label faceI,
            const label fp,
            const label v0,
            const label v1
        ) const;

        //- Gets labels of changed faces and enterVertices on faces.
        // Returns labels of faces changed and enterVertices on them.
        void faceToFace
        (
            const labelList& changedFaces,
            const labelList& enterVerts,

            labelList& nbrFaces,
            labelList& nbrEnterVerts
        );


        //- Disallow default bitwise copy construct
        walkPatch(const walkPatch&);

        //- Disallow default bitwise assignment
        void operator=(const walkPatch&);

public:

    ClassName("walkPatch");


    // Constructors

        //- Construct from components
        walkPatch
        (
            const primitivePatch& pp,
            const labelList& faceZone,  // Per face which zone it belongs to
            const bool reverse,         // Reverse walk
            const label faceI,          // Current face
            const label enterVertI,     // Vertex across which this face
                                        // is visited.
            boolList& visited
        );


    // Member Functions

        const dynamicLabelList& visitOrder() const
        {
            return visitOrder_;
        }

        const dynamicLabelList& indexInFace() const
        {
            return indexInFace_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
