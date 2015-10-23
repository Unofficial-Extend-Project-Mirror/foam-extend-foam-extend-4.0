/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    tetPolyBoundaryMesh

Description

SourceFiles
    tetPolyBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef tetPolyBoundaryMesh_H
#define tetPolyBoundaryMesh_H

#include "tetPolyPatchList.H"
#include "faceListFwd.H"
#include "lduInterfacePtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class tetPolyMesh;
class polyBoundaryMesh;
class globalTetPolyPatch;

/*---------------------------------------------------------------------------*\
                 Class tetPolyBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class tetPolyBoundaryMesh
:
    public tetPolyPatchList
{
    // Private data

        //- Reference to mesh
        const tetPolyMesh& mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPolyBoundaryMesh(const tetPolyBoundaryMesh&);

        //- Disallow default bitwise assignment
        void operator=(const tetPolyBoundaryMesh&);


public:

    // Constructors

        //- Construct from polyBoundaryMesh
        tetPolyBoundaryMesh
        (
            const tetPolyMesh&,
            const polyBoundaryMesh&
        );


    // Member functions

        // Access

        //- Return the mesh reference
            const tetPolyMesh& mesh() const
            {
                return mesh_;
            }

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            lduInterfacePtrsList interfaces() const;

            //- Return reference to shared processor point patch
            const globalTetPolyPatch& globalPatch() const;

            //- Return boundary faces per patch.  All are triangles
            faceListList boundaryTriFaces() const;


        // Edit

            // Update topology
            void updateMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
