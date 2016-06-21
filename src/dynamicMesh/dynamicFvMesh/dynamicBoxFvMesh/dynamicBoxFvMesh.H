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
    dynamicBoxFvMesh

Description
    Automatic simplified mesh motion for "box-in-mesh" cases.  Here,
    a direction of motion is defined, together with motion amplitude
    and frequency.  The domain is separated into three parts, where
    the middle part moves according to the prescribed motion law.
    Parts of the mesh before and after the obstacle are scaled.

SourceFiles
    dynamicBoxFvMesh.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef dynamicBoxFvMesh_H
#define dynamicBoxFvMesh_H

#include "dynamicFvMesh.H"
#include "dictionary.H"
#include "pointIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicBoxFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicBoxFvMesh
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

        //- Motion amplitude
        vector amplitude_;

        //- Motion frequency
        scalar frequency_;

        //- Stationary point location
        pointIOField stationaryPoints_;

        //- Vertex markup field
        scalarField motionMarkup_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicBoxFvMesh(const dynamicBoxFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicBoxFvMesh&);


public:

    //- Runtime type information
    TypeName("dynamicBoxFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicBoxFvMesh(const IOobject& io);


    // Destructor

        virtual ~dynamicBoxFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
