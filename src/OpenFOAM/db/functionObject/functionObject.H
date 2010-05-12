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
    Foam::functionObject

Description
    Abstract base-class for Time/database function objects.

SourceFiles
    functionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObject_H
#define functionObject_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;

/*---------------------------------------------------------------------------*\
                           Class functionObject Declaration
\*---------------------------------------------------------------------------*/

class functionObject
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        functionObject(const functionObject&);

        //- Disallow default bitwise assignment
        void operator=(const functionObject&);


public:

    //- Runtime type information
    virtual const word& type() const = 0;

    static int debug;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            functionObject,
            dictionary,
            (const word& name, const Time& t, const dictionary& dict),
            (name, t, dict)
        );


    // Constructors

        //- Construct null
        functionObject();

        //- Return clone
        autoPtr<functionObject> clone() const
        {
            notImplemented("functionObject::clone() const");
            return autoPtr<functionObject>(NULL);
        }

        //- Return a pointer to a new functionObject created on freestore
        //  from Istream
        class iNew
        {
            const Time& time_;

        public:

            iNew(const Time& t)
            :
                time_(t)
            {}

            autoPtr<functionObject> operator()
            (
                const word& name,
                Istream& is
            ) const;
        };


    // Selectors

        //- Select from Time and Istream
        static autoPtr<functionObject> New
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Destructor

        virtual ~functionObject();


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start() = 0;

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute() = 0;

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
