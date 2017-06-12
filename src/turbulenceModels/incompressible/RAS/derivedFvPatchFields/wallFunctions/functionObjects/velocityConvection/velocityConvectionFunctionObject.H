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
    velocityConvectionFunctionObject

Description
    FunctionObject that calculates the velocity convection given names for the
    velocity field and the flux field and registers it into the database. Useful
    for convection sensitised wall functions (see omega/nutCWT/MEWT wall
    functions).

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved

SourceFiles
    velocityConvectionFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef velocityConvectionFunctionObject_H
#define velocityConvectionFunctionObject_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "volFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class velocityConvectionFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class velocityConvectionFunctionObject
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Region name
        word regionName_;

        //- Name of the velocity field
        word UName_;

        //- Name of the flux field
        word phiName_;

        //- Convection field
        volVectorField convection_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        velocityConvectionFunctionObject
        (
            const velocityConvectionFunctionObject&
        );

        //- Disallow default bitwise assignment
        void operator=(const velocityConvectionFunctionObject&);


public:

    //- Runtime type information
    TypeName("velocityConvection");


    // Constructors

        //- Construct from components
        velocityConvectionFunctionObject
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute();

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
