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
    Foam::GeoMesh

Description
    Generic mesh wrapper used by volMesh, surfaceMesh, pointMesh etc.

\*---------------------------------------------------------------------------*/

#ifndef GeoMesh_H
#define GeoMesh_H

#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class GeoMesh Declaration
\*---------------------------------------------------------------------------*/

template<class MESH>
class GeoMesh
{

protected:

    // Protected data

        //- Reference to Mesh
        const MESH& mesh_;


public:

    // Public typedefs

        typedef MESH Mesh;
        typedef typename MESH::BoundaryMesh BoundaryMesh;


    // Constructors

        //- Construct from MESH
        explicit GeoMesh(const MESH& mesh)
        :
            mesh_(mesh)
        {}


    // Member Functions

        //- Return the object registry
        const objectRegistry& thisDb() const
        {
            return mesh_;
        }


    // Member Operators

        //- Return reference to polyMesh
        const MESH& operator()() const
        {
            return mesh_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
