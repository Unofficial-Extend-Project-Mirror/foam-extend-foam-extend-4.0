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
    eBoundaryMesh

SourceFiles
    eBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef eBoundaryMesh_H
#define eBoundaryMesh_H

#include "ePatchList.H"
#include "wordList.H"
#include "pointField.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class eMesh;

/*---------------------------------------------------------------------------*\
                      Class eBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class eBoundaryMesh
:
    public ePatchList,
    public regIOobject
{
    // private data

        //- Reference to mesh
        const eMesh& mesh_;

        //- Disallow construct as copy
        eBoundaryMesh(const eBoundaryMesh&);

        //- Disallow assignment
        void operator=(const eBoundaryMesh&);

public:

    //- Runtime type information
    TypeName("eBoundaryMesh");


    // Constructors

        //- Construct from dictionary
        eBoundaryMesh
        (
            const IOobject& io,
            const eMesh& mesh
        );

        //- Construct given size
        eBoundaryMesh
        (
            const IOobject& io,
            const eMesh& mesh,
            const label size
        );


    // Destructor - default


    // Member functions

        // Access

            //- Return the mesh reference
            const eMesh& mesh() const;

            //- Return a list of patch types
            wordList types() const;

            //- Return a list of patch names
            wordList names() const;

            //- Return a list of patch sizes
            labelList patchSizes() const;

            //- Return a list of patch starts
            labelList patchStarts() const;

            //- Return patch index for a given edge label
            label whichPatch(const label edgeIndex) const;

            //- Find patch index given a name
            label findPatchID(const word& patchName) const;

            //- writeData member function required by regIOobject
            bool writeData(Ostream&) const;

        // Input / Output

            //- Read from input stream
            void readFromInputStream();

    // Ostream operator

        friend Ostream& operator<<(Ostream&, const eBoundaryMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
