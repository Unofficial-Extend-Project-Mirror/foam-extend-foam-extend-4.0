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
    Foam::OutputFilterFunctionObject

Description
    FunctionObject wrapper around OutputFilter to allow them to be created
    via the functions list within controlDict.

SourceFiles
    OutputFilterFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef OutputFilterFunctionObject_H
#define OutputFilterFunctionObject_H

#include "functionObject.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class OutputFilterFunctionObject Declaration
\*---------------------------------------------------------------------------*/

template<class OutputFilter>
class OutputFilterFunctionObject
:
    public functionObject
{
    // Private data

        word name_;
        const Time& time_;
        dictionary dict_;
        word regionName_;
        word dictName_;

        //- the execution interval (in time steps)
        //  a value <= 1 means execute at every time step
        label interval_;

        //- Switch for the execution of the functionObjects
        bool execution_;

        autoPtr<OutputFilter> ptr_;


    // Private Member Functions

        //- Read relevant dictionary entries
        void readDict();

        //- Disallow default bitwise copy construct
        OutputFilterFunctionObject(const OutputFilterFunctionObject&);

        //- Disallow default bitwise assignment
        void operator=(const OutputFilterFunctionObject&);


public:

    //- Runtime type information
    TypeName(OutputFilter::typeName_());


    // Constructors

        //- Construct from components
        OutputFilterFunctionObject
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- Return name
        virtual const word& name() const
        {
            return name_;
        }

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute();

        //- Switch the function object on
        virtual void on();

        //- Switch the function object off
        virtual void off();

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "OutputFilterFunctionObject.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
