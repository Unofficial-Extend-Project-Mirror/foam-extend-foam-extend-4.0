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
    Foam::JobInfo

Description
    Helper class for recording information about run/finished jobs.

    Writes the following files:
      - $FOAM_JOB_DIR/runningJobs
      - $FOAM_JOB_DIR/finishedJobs

SourceFiles
    JobInfo.C

\*---------------------------------------------------------------------------*/

#ifndef JobInfo_H
#define JobInfo_H

#include "dictionary.H"
#include "fileName.H"
#include "cpuTime.H"
#include "infoSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class JobInfo Declaration
\*---------------------------------------------------------------------------*/

class JobInfo
:
    public dictionary
{
    // Private data
        fileName runningJobPath_;
        fileName finishedJobPath_;
        cpuTime  cpuTime_;

    // Private member functions

        bool write(Ostream&) const;
        void end(const word& terminationType);


public:

    static bool constructed;
    static debug::infoSwitch writeJobInfo;

    // Constructors

        //- Construct null
        JobInfo();


    // Destructor

        ~JobInfo();


    // Member Functions

        // Write

            void write() const;

            void end();

            void exit();

            void abort();

            void signalEnd() const;
};


extern JobInfo jobInfo;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
