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
    Foam::subCycleTime

Description
    A class for managing sub-cycling times

SourceFiles
    subCycleTime.C

\*---------------------------------------------------------------------------*/

#ifndef subCycleTime_H
#define subCycleTime_H

#include "objectRegistry.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class subCycleTime Declaration
\*---------------------------------------------------------------------------*/

class subCycleTime
{
    // Private data

        Time& time_;

        label nSubCycles_;
        label subCycleIndex_;

public:

    // Constructors

        //- Construct from original time and number of sub-cycles
        subCycleTime(Time&, const label nSubCycles);


    // Destructor

        ~subCycleTime();


    // Member functions

        //- Return true if the number of sub-cycles has been reached
        bool end() const;

        //- End the sub-cycling and reset the time-state
        void endSubCycle();


    // Member operators

        //- Prefix increment
        subCycleTime& operator++();

        //- Postfix increment
        subCycleTime& operator++(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
