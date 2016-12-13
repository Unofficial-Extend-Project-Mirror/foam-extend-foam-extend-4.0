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
    Foam::MeshObject

Description
    Templated abstract base-class for dynamic mesh objects used to automate
    their allocation to the mesh database and the mesh-modifier event-loop.

SourceFiles
    MeshObject.C

\*---------------------------------------------------------------------------*/

#ifndef MeshObject_H
#define MeshObject_H

#include "meshObjectBase.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class MeshObject Declaration
\*---------------------------------------------------------------------------*/


template<class Mesh, class Type>
class MeshObject
:
    public meshObjectBase,
    public regIOobject
{
    // Private data

        // Reference to Mesh
        const Mesh& mesh_;


public:

    // Constructors

        explicit MeshObject(const Mesh& mesh);

        static const Type& New(const Mesh& mesh);

        template<class Data1>
        static const Type& New
        (
            const Mesh& mesh,
            const Data1& d
        );

        template<class Data1, class Data2>
        static const Type& New
        (
            const Mesh& mesh,
            const Data1&,
            const Data2&
        );

        template<class Data1, class Data2, class Data3>
        static const Type& New
        (
            const Mesh& mesh,
            const Data1&,
            const Data2&,
            const Data3&
        );

        template<class Data1, class Data2, class Data3, class Data4>
        static const Type& New
        (
            const Mesh& mesh,
            const Data1&,
            const Data2&,
            const Data3&,
            const Data4&
        );


    // Destructors

        virtual ~MeshObject();

        static bool Delete(const Mesh& mesh);


    // Member Functions

        //- Delete object from database
        virtual bool deleteObject() const;

        //- Return mesh reference
        const Mesh& mesh() const
        {
            return mesh_;
        }

        virtual bool writeData(Foam::Ostream&) const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MeshObject.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
