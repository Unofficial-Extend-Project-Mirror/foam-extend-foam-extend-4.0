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
    dynamicBodyFvMesh

Description
    Automatic motion of the mesh around a moving body. A direction,
    amplitude and frequency of translational motion and origin, axis,
    amplitude and frequency of rotational motion must be specified.

SourceFiles
    dynamicBodyFvMesh.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef dynamicBodyFvMesh_H
#define dynamicBodyFvMesh_H

#include "dynamicFvMesh.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                   Class dynamicBodyFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicBodyFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Moving mesh dictionary
        dictionary dynamicMeshCoeffs_;

        //- Moving mesh solver
        autoPtr<motionSolver> motionPtr_;

        //- Name of the patch which encloses the body
        word bodyPatchName_;

        //- Body patch ID
        label bodyPatchID_;

        //- Translational motion direction
        vector translationDirection_;

        //- Translational motion amplitude
        scalar translationAmplitude_;

        //- Translational motion frequency
        scalar translationFrequency_;

        //- Initial rotational motion origin
        vector initialRotationOrigin_;

        //- Rotational motion axis
        vector rotationAxis_;

        //- Rotational motion amplitude
        scalar rotationAmplitude_;

        //- Rotational motion frequency
        scalar rotationFrequency_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicBodyFvMesh(const dynamicBodyFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicBodyFvMesh&);


public:

    //- Runtime type information
    TypeName("dynamicBodyFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicBodyFvMesh(const IOobject& io);


    // Destructor

        virtual ~dynamicBodyFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
