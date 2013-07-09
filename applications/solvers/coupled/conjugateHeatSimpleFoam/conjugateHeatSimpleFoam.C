/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 H. Jasak & H. Rusche
     \\/     M anipulation  | All rights reserved
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
    chipAirCoolingSimpleFoam

Description
    Steady-State solver for buoyancy-driven turbulent flow of incompressible
    Newtonian fluids with conjugate heat transfer, complex heat conduction
    and radiation

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "coupledFvMatrices.H"
#include "regionCouplePolyPatch.H"
#include "radiationModel.H"
#include "thermalModel.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createFluidMesh.H"
#   include "createSolidMesh.H"
#   include "readGravitationalAcceleration.H"
#   include "createFields.H"
#   include "createSolidFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSIMPLEControls.H"
#       include "initConvergenceCheck.H"

        // Detach patches
#       include "detachPatches.H"

        p_rgh.storePrevIter();

#       include "UEqn.H"
#       include "pEqn.H"

        // Update turbulent quantities
        turbulence->correct();

        radiation->correct();

        // Update thermal conductivity in the fluid
        kappaEff = rho*Cp*(turbulence->nu()/Pr + turbulence->nut()/Prt);

        // Update thermal conductivity in the solid
        solidThermo.correct();
        ksolid = solidThermo.k();

        // Coupled patches
#       include "attachPatches.H"

        kappaEff.correctBoundaryConditions();
        ksolid.correctBoundaryConditions();

        // Interpolate to the faces and add thermal resistance
        surfaceScalarField ksolidf = fvc::interpolate(ksolid);
        solidThermo.modifyResistance(ksolidf);

#       include "solveEnergy.H"

        // Update density according to Boussinesq approximation
        rhok = 1.0 - beta*(T - TRef);

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
