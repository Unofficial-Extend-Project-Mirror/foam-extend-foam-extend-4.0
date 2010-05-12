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
    Foam::functionObjectList

Description
    List of function objects with execute function which is called for
    each object.

SourceFiles
    functionObjectList.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjectList_H
#define functionObjectList_H

#include "functionObject.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class functionObjectList Declaration
\*---------------------------------------------------------------------------*/

class functionObjectList
:
    public HashPtrTable<functionObject>
{
    // Private data

        const Time& time_;

        //- Dictionary containing the list of function object specifications
        const dictionary& foDict_;

        //- Switch for the execultion of the functionObjects
        bool execution_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        functionObjectList(const functionObjectList&);

        //- Disallow default bitwise assignment
        void operator=(const functionObjectList&);


public:

    // Constructors

        //- Construct from Time and the execution setting
        //  The functionObject specifications are read from the controlDict
        functionObjectList
        (
            const Time& t,
            const bool execution = true
        );


        //- Construct from Time, functionObject dictionary and the execution
        //  setting
        functionObjectList
        (
            const Time& t,
            const dictionary& foDict,
            const bool execution = true
        );


    // Destructor

        virtual ~functionObjectList();


    // Member Functions

        //- Start is called at the start of the time-loop
        virtual bool start();

        //- Execute is called at each ++ or += of the time-loop
        virtual bool execute();

        //- Switch the function objects on
        virtual void on();

        //- Switch the function objects off
        virtual void off();

        //- Read and set the function objects if their data has changed
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
