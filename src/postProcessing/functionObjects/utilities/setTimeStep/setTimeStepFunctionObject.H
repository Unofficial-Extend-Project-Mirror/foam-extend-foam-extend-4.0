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
    Foam::setTimeStepFunctionObject

Group
    grpUtilitiesFunctionObjects

Description
    Overrides the timeStep. Can only be used with
    solvers with adjustTimeStep control (e.g. pimpleFoam). Makes no attempt
    to cooperate with other timeStep 'controllers' (maxCo, other
    functionObjects). Supports 'enabled' flag but none of othe other ones
    'timeStart', 'timeEnd', 'outputControl' etc.

SourceFiles
    setTimeStepFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef setTimeStepFunctionObject_H
#define setTimeStepFunctionObject_H

#include "functionObject.H"
#include "dictionary.H"
#include "Function1.H"
#include "pointFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class setTimeStepFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class setTimeStepFunctionObject
:
    public functionObject
{
    // Private data

        //- Reference to the time database
        const Time& time_;


        // Optional user inputs

            //- Switch for the execution - defaults to 'yes/on'
            bool enabled_;

            //- Time step
            autoPtr<Function1<scalar> > timeStepPtr_;


        //- Disallow default bitwise copy construct
        setTimeStepFunctionObject(const setTimeStepFunctionObject&);

        //- Disallow default bitwise assignment
        void operator=(const setTimeStepFunctionObject&);


public:
    //- Runtime type information
    TypeName("setTimeStep");

    // Constructors

        //- Construct from components
        setTimeStepFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Member Functions

        // Access

            //- Return time database
            virtual const Time& time() const
            {
                return time_;
            }

            //- Return the enabled flag
            virtual bool enabled() const
            {
                return enabled_;
            }


        // Function object control

            //- Switch the function object on
            virtual void on();

            //- Switch the function object off
            virtual void off();


            //- Called at the start of the time-loop
            virtual bool start();

            //- Called at each ++ or += of the time-loop
            virtual bool execute(const bool forceWrite);

            //- Called when Time::run() determines that the time-loop exits
            virtual bool end();

            //- Called when time was set at the end of the Time::operator++
            virtual bool timeSet();

            //- Called at the end of Time::adjustDeltaT() if adjustTime is true
            virtual bool adjustTimeStep();

            //- Read and set the function object if its data have changed
            virtual bool read(const dictionary&);

            //- Update for changes of mesh
            virtual void updateMesh(const mapPolyMesh& mpm);

            //- Update for changes of mesh
            virtual void movePoints(const pointField& points);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
