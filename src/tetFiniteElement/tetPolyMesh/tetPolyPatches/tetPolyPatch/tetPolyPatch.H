/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    tetPolyPatch

Description

SourceFiles
    tetPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef tetPolyPatch_H
#define tetPolyPatch_H

#include "labelList.H"
#include "vectorField.H"
#include "triFaceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class tetPolyMesh;
class tetPolyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                    Class tetPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class tetPolyPatch
{
    // Private data

        //- Reference to boundary mesh
        const tetPolyBoundaryMesh& boundaryMesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPolyPatch(const tetPolyPatch&);

        //- Disallow default bitwise assignment
        void operator=(const tetPolyPatch&);


public:

    typedef tetPolyMesh Mesh;


    //- Runtime type information
    TypeName("basePatch");


    // Constructor

        //- Construct from polyPatch
        tetPolyPatch
        (
            const tetPolyBoundaryMesh& bm
        )
        :
            boundaryMesh_(bm)
        {}


    // Destructor

        virtual ~tetPolyPatch()
        {}


    // Member Functions

        //- Return name
        virtual const word& name() const = 0;

        //- Return size
        virtual label size() const = 0;

        //- Return number of faces
        virtual label nFaces() const = 0;

        //- Return the index of this patch in the tetPolyBoundaryMesh
        virtual label index() const = 0;


        //- Return boundaryMesh reference
        const tetPolyBoundaryMesh& boundaryMesh() const
        {
            return boundaryMesh_;
        }

        //- Return mesh points
        virtual const labelList& meshPoints() const = 0;

        //- Return mesh points
        virtual const vectorField& localPoints() const = 0;

        //- Return  point normals
        virtual const vectorField& pointNormals() const = 0;

        //- Face decomposition into triangle
        //  Used for implementation of FEM boundary conditions
        virtual triFaceList faceTriangles(const label faceID) const = 0;

        //- Return complete list of faces. All are triangles
        virtual faceList triFaces() const = 0;


        // Update topology
        virtual void updateMesh()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
