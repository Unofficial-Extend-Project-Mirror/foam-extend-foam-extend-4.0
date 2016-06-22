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
    Foam::cpuTime

Description
    Starts timing CPU usage and return elapsed time from start.

SeeAlso
    clockTime

SourceFiles
    cpuTime.C

\*---------------------------------------------------------------------------*/

#ifndef cpuTime_H
#define cpuTime_H

#include <time.h>
#include <sys/times.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cpuTime Declaration
\*---------------------------------------------------------------------------*/

class cpuTime
{
    // Private data

        static long Hz_;

        struct tms startTime_;
        mutable struct tms lastTime_;
        mutable struct tms newTime_;

        static void getTime(struct tms& t);

        static double timeDifference
        (
            const struct tms& start,
            const struct tms& end
        );


public:

    // Constructors

        //- Construct from components
        cpuTime();


    // Member Functions

        // Access

            //- Returns CPU time from start of run
            double elapsedCpuTime() const;

            //- Returns CPU time from last call of cpuTimeIncrement()
            double cpuTimeIncrement() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
