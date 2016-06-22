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
    Foam::solution

Description
    Selector class for relaxation factors, solver type and solution.

SourceFiles
    solution.C

\*---------------------------------------------------------------------------*/

#ifndef solution_H
#define solution_H

#include "IOdictionary.H"
#include "debugSwitch.H"
#include "lduSolverPerformance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class solution Declaration
\*---------------------------------------------------------------------------*/

class solution
:
    public IOdictionary
{
    // Private data

        //- Dictionary of temporary fields to cache
        dictionary cache_;

        //- Switch for the caching mechanism
        bool caching_;

        //- Dictionary of relaxation factors for all the fields
        dictionary fieldRelaxDict_;

        //- Dictionary of relaxation factors for all the equations
        dictionary eqnRelaxDict_;

        //- Optional default relaxation factor for all the fields
        scalar fieldRelaxDefault_;

        //- Optional default relaxation factor for all the equations
        scalar eqnRelaxDefault_;

        //- Dictionary of solver parameters for all the fields
        dictionary solvers_;

        //- Dictionary of solver performance.  Used in convergence check
        mutable dictionary solverPerformance_;

        //- Previously used time-index, used for reset between iterations
        mutable label prevTimeIndex_;


    // Private Member Functions

        //- Read settings from the dictionary
        void read(const dictionary&);

        //- Disallow default bitwise copy construct and assignment
        solution(const solution&);
        void operator=(const solution&);


public:

    //- Update from older solver controls syntax
    //  Usually verbose, since we want to know about the changes
    //  Returns the number of settings changed
    static label upgradeSolverDict(dictionary& dict, const bool verbose=true);

    //- Debug switch
    static debug::debugSwitch debug;


    // Constructors

        //- Construct for given objectRegistry and dictionary
        solution
        (
            const objectRegistry& obr,
            const fileName& dictName
        );


    // Member Functions

        // Access

            //- Return true if the given field should be cached
            bool cache(const word& name) const;

            //- Helper for printing cache message
            template<class FieldType>
            static void cachePrintMessage
            (
                const char* message,
                const word& name,
                const FieldType& vf
            );

            //- Return true if the relaxation factor is given for the field
            bool relaxField(const word& name) const;

            //- Return true if the relaxation factor is given for the equation
            bool relaxEquation(const word& name) const;

            //- Return the relaxation factor for the given field
            scalar fieldRelaxationFactor(const word& name) const;

            //- Return the relaxation factor for the given eqation
            scalar equationRelaxationFactor(const word& name) const;

            //- Return the selected sub-dictionary of solvers if the "select"
            //  keyword is given, otherwise return the complete dictionary
            const dictionary& solutionDict() const;

            //- Return the solver controls dictionary for the given field
            const dictionary& solverDict(const word& name) const;

            //- Return the solver controls dictionary for the given field
            const dictionary& solver(const word& name) const;

            //- Return the dictionary of solver performance data
            //  which includes initial and final residuals for convergence
            //  checking
            dictionary& solverPerformanceDict() const;

            //- Add/set the solverPerformance entry for the named field
            void setSolverPerformance
            (
                const word& name,
                const lduSolverPerformance&
            ) const;

            //- Add/set the solverPerformance entry, using its fieldName
            void setSolverPerformance
            (
                const lduSolverPerformance&
            ) const;


        // Edit

            //- Return access to field relaxation factors dictionary
            dictionary& fieldRelaxationFactors()
            {
                return fieldRelaxDict_;
            }

            //- Return access to equation relaxation factors dictionary
            dictionary& equationRelaxationFactors()
            {
                return eqnRelaxDict_;
            }

            //- Return access to solvers dictionary
            dictionary& solvers()
            {
                return solvers_;
            }


        // Read

            //- Read the solution dictionary
            bool read();


        // Write

            //- WriteData function required for regIOobject write operation
            virtual bool writeData(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "solutionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
