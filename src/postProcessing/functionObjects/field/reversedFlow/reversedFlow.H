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
    reversedFlow

Description
    FunctionObject detects reversed flow on an (outlet) patch

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    reversedFlow.C

\*---------------------------------------------------------------------------*/

#ifndef reversedFlow_H
#define reversedFlow_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class reversedFlow Declaration
\*---------------------------------------------------------------------------*/

class reversedFlow
:
    public functionObject
{
    // Private data

        //- Reference to main object registry
        const Time& time_;

        //- Region name
        word regionName_;

        //- Patch name
        word patchName_;

        //- Flux field name
        word fluxFieldName_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        reversedFlow
        (
            const reversedFlow&
        );

        //- Disallow default bitwise assignment
        void operator=(const reversedFlow&);


public:

    //- Runtime type information
    TypeName("reversedFlow");


    // Constructors

        //- Construct from components
        reversedFlow
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    //- Destructor
    virtual ~reversedFlow()
    {}


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute(const bool forceWrite);

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for mesh point-motion
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
