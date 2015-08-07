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

Application
    sixDOFsolver

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.

Description
    Six degrees of freedom solver for multiple bodies

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "ODESolver.H"
#include "sixDOFbodies.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    sixDOFbodies structure(runTime);
    OFstream of(runTime.path()/"motion.dat");

    Info<< "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info << "Time = " << runTime.value() << endl;

        structure.solve();

        of << runTime.value() << tab;

        forAll (structure.names(), bodyI)
        {
            Info<< nl << "Average velocity of " << structure.names()[bodyI]
                << " in time step = "
                << structure()[bodyI].Uaverage().value() << nl
                << "Current velocity in time instant = "
                << structure()[bodyI].U().value() << nl
                << "Average omega of " << structure.names()[bodyI]
                << " in time step = "
                << structure()[bodyI].omegaAverage().value() << nl
                << "Current omega in time instant = "
                << structure()[bodyI].omega().value()  << nl
                << "Average omegaAbs in time step = "
                << structure()[bodyI].omegaAverageAbsolute().value() << nl
                << endl;

            of  << structure()[bodyI].X().value().x() << tab;
        }

        of << endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
