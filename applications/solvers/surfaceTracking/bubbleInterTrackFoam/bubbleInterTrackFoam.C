/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Z. Tukovic and H. Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    interfaceTrackinFoam

Description
    Incompressible laminar CFD code for simulation of a single bubble rising 
    in a stil liquid. Interface between fluid phases is tracked using moving
    mesh.

\*---------------------------------------------------------------------------*/

#include "Ostream.H"
#include "OFstream.H"

#include "fvCFD.H"

#include "inletOutletFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

#include "motionSolver.H"
#include "freeSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"

#   include "createMeshNoClear.H"

#   include "createFields.H"

#   include "initContinuityErrs.H"

#   include "createBubble.H"

#   include "createPath.H"

#   include "createForces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info << "Time = " << runTime.value() << endl;

#       include "readBubbleSIMPLEControls.H"
#       include "CourantNo.H"

        interface.updateDisplacementDirections();

        interface.moveMeshPointsForOldFreeSurfDisplacement();

        interface.smooth();

#       include "updateMovingReferenceFrame.H"

        // --- SIMPLE loop
        for (label timeCorr=0; timeCorr<=nTimeCorr; timeCorr++)
        {
            Info << endl;

            p.storePrevIter();        

            interface.correctBoundaryConditions();

            tmp<fvVectorMatrix> UEqn
            (
                fvm::ddt(rho, U)
              + rho*aF
              + fvm::div(phiNet, U)
              - fvm::laplacian(mu, U)             
            );

            UEqn().relax();

            solve(UEqn() == -fvc::grad(p));

            volScalarField AU = UEqn().A();

            U = UEqn().H()/AU;
            U.correctBoundaryConditions();

            UEqn.clear();

            phi = (fvc::interpolate(U) & mesh.Sf());

#           include "spacePatchFlowRateScaling.H"

#           include "freeSurfacePatchFlowRateScaling.H"

            // Non-orthogonal pressure corrector loop
            for (label nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(1.0/AU, p) == fvc::div(phi)
                );

                pEqn.setReference(0, 0.0);
                pEqn.solve();

                if (nonOrth == nNonOrthCorr)
                {
                    phi -= pEqn.flux();
                }
            }

#           include "continuityErrs.H"

            // Explicitly relax pressure
            p.relax();

            // Momentum corrector
            U -= fvc::grad(p)/AU;
            U.correctBoundaryConditions();

            // Move mesh
            interface.movePoints();

            // Update motion fluxes
            phiNet = fvc::interpolate(rho)*(phi - fvc::meshPhi(rho, U));

#           include "freeSurfaceContinuityErrs.H"
        }

#       include "bubbleVolumeReduction.H"

        runTime.write();

#       include "writeData.H"

#       include "writePath.H"

#       include "writeForces.H"

        Info << "ExecutionTime = "
            << scalar(runTime.elapsedCpuTime())
            << " s\n" << endl << endl;
    }

    Info << "End\n" << endl;

    return(0);
}

// ************************************************************************* //
