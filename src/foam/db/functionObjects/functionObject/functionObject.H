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
    Foam::functionObject

Description
    Abstract base-class for Time/database function objects.

See Also
    Foam::OutputFilterFunctionObject

SourceFiles
    functionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObject_H
#define functionObject_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "debugSwitch.H"
#include "pointFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class functionObject Declaration
\*---------------------------------------------------------------------------*/

class functionObject
{
    // Private data

        //- Name
        const word name_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        functionObject(const functionObject&);

        //- Disallow default bitwise assignment
        void operator=(const functionObject&);


public:

    //- Runtime type information
    virtual const word& type() const = 0;

    static debug::debugSwitch debug;


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

        //- Construct from components
        functionObject(const word& name);

        //- Return clone
        autoPtr<functionObject> clone() const
        {
            notImplemented("functionObject::clone() const");
            return autoPtr<functionObject>(nullptr);
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

        //- Select from dictionary, based on its "type" entry
        static autoPtr<functionObject> New
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    //- Destructor

        virtual ~functionObject();


    // Member Functions

        //- Name
        virtual const word& name() const;

        //- Called at the start of the time-loop
        virtual bool start() = 0;

        //- Called at each ++ or += of the time-loop. forceWrite overrides the
        //  outputControl behaviour.
        virtual bool execute(const bool forceWrite) = 0;

        //- Called when Time::run() determines that the time-loop exits.
        //  By default it simply calls execute().
        virtual bool end();

        //- Called when time was set at the end of the Time::operator++
        virtual bool timeSet();

        //- Called at the end of Time::adjustDeltaT() if adjustTime is true
        virtual bool adjustTimeStep();

        //- Read and set the function object if its data have changed
        virtual bool read(const dictionary&) = 0;

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh& mpm) = 0;

        //- Update for changes of mesh
        virtual void movePoints(const pointField& points) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
