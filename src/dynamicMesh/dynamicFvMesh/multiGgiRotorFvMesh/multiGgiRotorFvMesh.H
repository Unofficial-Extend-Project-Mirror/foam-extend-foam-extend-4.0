/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    multiGgiRotorFvMesh

Description
    Mixer mesh with multiple rotors.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    multiGgiRotorFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef multiGgiRotorFvMesh_H
#define multiGgiRotorFvMesh_H

#include "dynamicFvMesh.H"
#include "ggiRotor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class multiGgiRotorFvMesh Declaration
\*---------------------------------------------------------------------------*/

class multiGgiRotorFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary dict_;

        //- Mixer rotors
        PtrList<ggiRotor> rotors_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiGgiRotorFvMesh(const multiGgiRotorFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const multiGgiRotorFvMesh&);


public:

    //- Runtime type information
    TypeName("multiGgiRotorFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit multiGgiRotorFvMesh(const IOobject& io);


    //- Destructor
    virtual ~multiGgiRotorFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
