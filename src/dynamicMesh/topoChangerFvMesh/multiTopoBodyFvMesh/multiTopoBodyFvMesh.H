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
    multiTopoBodyFvMesh

Description
    Mesh with multiple moving bodies using topological layering.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    multiTopoBodyFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef multiTopoBodyFvMesh_H
#define multiTopoBodyFvMesh_H

#include "topoChangerFvMesh.H"
#include "topoBody.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class multiTopoBodyFvMesh Declaration
\*---------------------------------------------------------------------------*/

class multiTopoBodyFvMesh
:
    public topoChangerFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary dict_;

        //- Moving bodies
        PtrList<topoBody> bodies_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiTopoBodyFvMesh(const multiTopoBodyFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const multiTopoBodyFvMesh&);


        //- Add mixer rotor zones and modifiers
        void addZonesAndModifiers();

        //- Sum up point motion from all topo bodies
        tmp<vectorField> pointMotion() const;


public:

    //- Runtime type information
    TypeName("multiTopoBodyFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit multiTopoBodyFvMesh(const IOobject& io);


    // Destructor

        virtual ~multiTopoBodyFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
