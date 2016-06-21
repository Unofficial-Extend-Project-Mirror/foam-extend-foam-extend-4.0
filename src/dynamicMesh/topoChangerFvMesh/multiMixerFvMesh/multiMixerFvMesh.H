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
    multiMixerFvMesh

Description
    Mixer mesh with multiple rotors.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    multiMixerFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef multiMixerFvMesh_H
#define multiMixerFvMesh_H

#include "topoChangerFvMesh.H"
#include "mixerRotor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class multiMixerFvMesh Declaration
\*---------------------------------------------------------------------------*/

class multiMixerFvMesh
:
    public topoChangerFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary dict_;

        //- Mixer rotors
        PtrList<mixerRotor> rotors_;

        //- Attach-detach action
        Switch attachDetach_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiMixerFvMesh(const multiMixerFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const multiMixerFvMesh&);


        //- Add mixer rotor zones and modifiers
        void addZonesAndModifiers();

        //- Return true if any of the rotors uses a topological change
        bool useTopoSliding() const;

        //- Check that rotors using a topological change are in sync
        void checkRotors() const;

        //- Return true if sliding interface are attached
        bool attached() const;


public:

    //- Runtime type information
    TypeName("multiMixerFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit multiMixerFvMesh(const IOobject& io);


    // Destructor

        virtual ~multiMixerFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
