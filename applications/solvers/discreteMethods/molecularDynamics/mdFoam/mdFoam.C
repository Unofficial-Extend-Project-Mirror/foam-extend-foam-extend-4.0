/*---------------------------------------------------------------------------*\
  =========                   |
  \\      /   F ield          | foam-extend: Open Source CFD
   \\    /    O peration      |
    \\  /     A nd            | For copyright notice see file Copyright
     \\/      M anipulation   |
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
    mdFoam

Description
    Molecular dynamics solver for fluid dynamics

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "md.H"

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    potential pot(mesh);

    moleculeCloud molecules(mesh, pot);

#   include "temperatureAndPressureVariables.H"

    label nAveragingSteps = 0;

    Info << "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {

        nAveragingSteps++;

        Info << "Time = " << runTime.timeName() << endl;

        molecules.evolve();

#       include "meanMomentumEnergyAndNMols.H"

#       include "temperatureAndPressure.H"

        runTime.write();

        if (runTime.outputTime())
        {
            nAveragingSteps = 0;
        }

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info << "End\n" << endl;

    return 0;
}
