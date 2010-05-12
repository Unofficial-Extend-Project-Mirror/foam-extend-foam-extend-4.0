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
    sonicFoam

Description
    Transient solver for trans-sonic/supersonic, turbulent flow of a
    compressible gas.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicThermo.H"
#include "compressible/RASModel/RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readPISOControls.H"
#       include "compressibleCourantNo.H"

#       include "rhoEqn.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(rho, U)
          + fvm::div(phi, U)
          + turbulence->divDevRhoReff(U)
        );

        solve(UEqn == -fvc::grad(p));

#       include "hEqn.H"


        // --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {
            volScalarField rUA = 1.0/UEqn.A();
            U = rUA*UEqn.H();

            surfaceScalarField phid
            (
                "phid",
                fvc::interpolate(thermo->psi())
               *(
                   (fvc::interpolate(U) & mesh.Sf())
                 + fvc::ddtPhiCorr(rUA, rho, U, phi)
               )
            );

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::ddt(psi, p)
                  + fvm::div(phid, p, "div(phid,p)")
                  - fvm::laplacian(rho*rUA, p)
                );

                pEqn.solve();

                if (nonOrth == nNonOrthCorr)
                {
                    phi = pEqn.flux();
                }
            }

#           include "compressibleContinuityErrs.H"

            U -= rUA*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        DpDt =
            fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

        turbulence->correct();

        rho = psi*p;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
