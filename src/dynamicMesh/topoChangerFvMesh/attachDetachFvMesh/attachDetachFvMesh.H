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
    attachDetachFvMesh

Description

SourceFiles
    attachDetachFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef attachDetachFvMesh_H
#define attachDetachFvMesh_H

#include "topoChangerFvMesh.H"
#include "motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class attachDetachFvMesh Declaration
\*---------------------------------------------------------------------------*/

class attachDetachFvMesh
:
    public topoChangerFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary motionDict_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        attachDetachFvMesh(const attachDetachFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const attachDetachFvMesh&);


        //- Add attach/detach zones and modifiers
        void addZonesAndModifiers();


public:

    //- Runtime type information
    TypeName("attachDetachFvMesh");


    // Constructors

        //- Construct from database
        explicit attachDetachFvMesh(const IOobject& io);


    // Destructor

        virtual ~attachDetachFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
