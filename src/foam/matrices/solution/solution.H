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

        //- Dictionary of relaxation factors for all the fields
        dictionary relaxationFactors_;

        //- Optional default relaxation factor for all the fields
        scalar defaultRelaxationFactor_;

        //- Dictionary of solver parameters for all the fields
        dictionary solvers_;


    // Private Member Functions

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
        solution(const objectRegistry& obr, const fileName& dictName);


    // Member Functions

        // Access

            //- Return true if the relaxation factor is given for the field
            bool relax(const word& name) const;

            //- Return the relaxation factor for the given field
            scalar relaxationFactor(const word& name) const;

            //- Return the solver controls dictionary for the given field
            const dictionary& solverDict(const word& name) const;

            //- Return the solver controls dictionary for the given field
            const dictionary& solver(const word& name) const;


        // Edit

            //- Return access to relaxation factors dictionary
            dictionary& relaxationFactors()
            {
                return relaxationFactors_;
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

#endif

// ************************************************************************* //
