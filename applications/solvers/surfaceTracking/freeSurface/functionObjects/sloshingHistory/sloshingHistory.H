/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    sloshingHistory

Description
    FunctionObject reports min and max of a volScalarField

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    sloshingHistory.C

\*---------------------------------------------------------------------------*/

#ifndef sloshingHistory_H
#define sloshingHistory_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class sloshingHistory Declaration
\*---------------------------------------------------------------------------*/

class sloshingHistory
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

        //- Output file stream
        OFstream* historyFilePtr_;

        //- Free-surface patch ID
        label freeSurfacePatchID_;

        //- Free-surface point at the left side wall
        label leftPointID_;

        //- Free-surface point at the right side wall
        label rightPointID_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        sloshingHistory
        (
            const sloshingHistory&
        );

        //- Disallow default bitwise assignment
        void operator=(const sloshingHistory&);


public:

    //- Runtime type information
    TypeName("sloshingHistory");


    // Constructors

        //- Construct from components
        sloshingHistory
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
