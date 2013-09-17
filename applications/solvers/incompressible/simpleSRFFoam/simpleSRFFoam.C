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
    simpleSRFFoam

Description
    Steady-state solver for incompressible, turbulent flow of non-Newtonian
    fluids with single rotating frame.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"
#include "SRFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

    //mesh.clearPrimitives();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSIMPLEControls.H"

        p.storePrevIter();

        // Pressure-velocity SIMPLE corrector
        {
            // Momentum predictor
            tmp<fvVectorMatrix> UrelEqn
            (
                fvm::div(phi, Urel)
              + turbulence->divDevReff(Urel)
              + SRF->Su()
            );

            UrelEqn().relax();

            solve(UrelEqn() == -fvc::grad(p));

            p.boundaryField().updateCoeffs();
            volScalarField AUrel = UrelEqn().A();
            Urel = UrelEqn().H()/AUrel;
            UrelEqn.clear();
            phi = fvc::interpolate(Urel) & mesh.Sf();
            adjustPhi(phi, Urel, p);

            // Non-orthogonal pressure corrector loop
            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(1.0/AUrel, p) == fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (nonOrth == nNonOrthCorr)
                {
                    phi -= pEqn.flux();
                }
            }

#           include "continuityErrs.H"

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Momentum corrector
            Urel -= fvc::grad(p)/AUrel;
            Urel.correctBoundaryConditions();
        }

        turbulence->correct();

        // Recalculate Uabs
        Uabs = Urel + SRF->U();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
