/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createEngineTime.H"
#   include "createEngineDynamicMesh.H"
#   include "readPIMPLEControls.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"
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
#       include "readPIMPLEControls.H"
#       include "checkTotalVolume.H"
#       include "readEngineTimeControls.H"
#       include "compressibleCourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "Crank angle = " << runTime.theta() << " CA-deg" << endl;

        // Make flux absolute
        phi += meshFlux;

        bool meshChanged = mesh.update();

#       include "volContinuity.H"

        mesh.setBoundaryVelocity(U);

        if (meshChanged)
        {
            thermo.correct();
            rho = thermo.rho();
            rho.correctBoundaryConditions();
        }

        meshFlux = fvc::interpolate(rho)*fvc::meshPhi(rho, U);

        phi = fvc::interpolate(rho)
            *((fvc::interpolate(U) & mesh.Sf()) - fvc::meshPhi(rho, U));

        DpDt = dpdt + fvc::div(phi/fvc::interpolate(rho), p)
            - fvc::div(phi/fvc::interpolate(rho) + fvc::meshPhi(U))*p;

        {
#           include "compressibleCourantNo.H"
        }

        // Pressure-velocity corrector
        int oCorr = 0;
        do
        {
#           include "rhoEqn.H"
#           include "UEqn.H"

            // --- PISO loop
            for (int corr = 1; corr <= nCorr; corr++)
            {
#               include "pEqn.H"
#               include "hEqn.H"
            }
        } while (++oCorr < nOuterCorr);

        turbulence->correct();

#       include "logSummary.H"

        rho = thermo.rho();

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
