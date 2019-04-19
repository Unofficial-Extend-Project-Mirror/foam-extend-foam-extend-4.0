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
    pressureGradientFunctionObject

Description
    FunctionObject that calculates the pressure gradient given a name for the
    pressure field and registers it into the database. Useful for pressure
    sensitised wall functions (see omega/nutCWT/MEWT wall functions).

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved

SourceFiles
    pressureGradientFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef pressureGradientFunctionObject_H
#define pressureGradientFunctionObject_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "volFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class pressureGradientFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class pressureGradientFunctionObject
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

        //- Name of the pressure field
        word pName_;

        //- Pressure gradient field
        volVectorField gradp_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        pressureGradientFunctionObject
        (
            const pressureGradientFunctionObject&
        );

        //- Disallow default bitwise assignment
        void operator=(const pressureGradientFunctionObject&);


public:

    //- Runtime type information
    TypeName("pressureGradient");


    // Constructors

        //- Construct from components
        pressureGradientFunctionObject
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute(const bool forceWrite);

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
