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
    movingBoxFvMesh

Description
    Automatic simplified mesh motion for "box-in-mesh" cases.  Here,
    a direction of motion is defined, together with motion amplitude
    and frequency.  The domain is separated into three parts, where
    the middle part moves according to the prescribed motion law.
    Parts of the mesh before and after the obstacle are scaled.

SourceFiles
    movingBoxFvMesh.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef movingBoxFvMesh_H
#define movingBoxFvMesh_H

#include "dynamicFvMesh.H"
#include "dictionary.H"
#include "pointIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class movingBoxFvMesh Declaration
\*---------------------------------------------------------------------------*/

class movingBoxFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Moving mesh dictionary
        dictionary movingMeshCoeffs_;

        //- Motion box and direction
        vector splitDirection_;

        //- Left edge of box
        vector leftEdge_;

        //- Right edge of box
        vector rightEdge_;

        //- Motion velocity
        vector velocity_;

        //- Stationary point location
        pointIOField stationaryPoints_;

        //- Vertex markup field
        scalarField motionMarkup_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        movingBoxFvMesh(const movingBoxFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const movingBoxFvMesh&);


public:

    //- Runtime type information
    TypeName("movingBoxFvMesh");


    // Constructors

        //- Construct from IOobject
        movingBoxFvMesh(const IOobject& io);


    // Destructor

        virtual ~movingBoxFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
