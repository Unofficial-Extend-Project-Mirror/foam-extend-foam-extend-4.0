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

Application
    sonicTurbDyMEngineFoam

Description
    Solver for compressible cold flow in internal combustion engines
    with mesh motion and topological changes.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "engineTime.H"
#include "dynamicFvMesh.H"
#include "engineTopoChangerMesh.H"
#include "basicPsiThermo.H"
#include "turbulenceModel.H"
#include "Switch.H"
#include "OFstream.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createEngineTime.H"
#   include "createEngineDynamicMesh.H"

    pimpleControl pimple(mesh);

#   include "createFields.H"
#   include "initContinuityErrs.H"
#   include "createTimeControls.H"
#   include "readEngineTimeControls.H"
#   include "CourantNo.H"
#   include "setInitialDeltaT.H"
#   include "startSummary.H"
#   include "createEngineOutput.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    thermo.correct();

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "checkTotalVolume.H"
#       include "readEngineTimeControls.H"
#       include "compressibleCourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "Crank angle = " << runTime.theta() << " CA-deg" << endl;

        // Make the fluxes absolute (using the ddt(rho, U) scheme)
        phi += fvc::interpolate(rho)*fvc::meshPhi(rho, U);

        bool meshChanged = mesh.update();

#       include "volContinuity.H"

        // Make the fluxes relative (using the ddt(rho, U) scheme)
        phi -= fvc::interpolate(rho)*fvc::meshPhi(rho, U);

        if (meshChanged)
        {
            thermo.correct();
            rho = thermo.rho();
            rho.correctBoundaryConditions();
        }

        if (meshChanged)
        {
#           include "compressibleCourantNo.H"
        }

        // Pressure-velocity corrector
        while (pimple.loop())
        {
#           include "rhoEqn.H"
#           include "hEqn.H"
#           include "UEqn.H"

            // --- PISO loop
            while (pimple.correct())
            {
#               include "pEqn.H"
            }

            turbulence->correct();
        }

#       include "logSummary.H"

        runTime.write();

#       include "infoDataOutput.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
