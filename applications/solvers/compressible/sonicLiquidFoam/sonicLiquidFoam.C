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
    sonicLiquidFoam

Description
    Transient solver for trans-sonic/supersonic, laminar flow of a
    compressible liquid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readThermodynamicProperties.H"
#   include "readTransportProperties.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readPIMPLEControls.H"
#       include "compressibleCourantNo.H"

#       include "rhoEqn.H"

        // --- PIMPLE loop
        label oCorr = 0;
        do
        {
            fvVectorMatrix UEqn
            (
                fvm::ddt(rho, U)
              + fvm::div(phi, U)
              - fvm::laplacian(mu, U)
            );

            solve(UEqn == -fvc::grad(p));

            // --- PISO loop
            for (int corr = 0; corr < nCorr; corr++)
            {
                volScalarField rAU("rAU", 1.0/UEqn.A());
                surfaceScalarField rhorAUf
                (
                    "rhorAUf",
                    fvc::interpolate(rho*rAU)
                );

                U = rAU*UEqn.H();

                surfaceScalarField phid
                (
                    "phid",
                    psi*
                    (
                        (fvc::interpolate(U) & mesh.Sf())
                      + fvc::ddtPhiCorr(rAU, rho, U, phi)
                    )
                );

                phi = (rhoO/psi)*phid;

                fvScalarMatrix pEqn
                (
                    fvm::ddt(psi, p)
                  + fvc::div(phi)
                  + fvm::div(phid, p)
                  - fvm::laplacian(rhorAUf, p)
                );

                pEqn.solve();

                phi += pEqn.flux();

#               include "rhoEqn.H"
#               include "compressibleContinuityErrs.H"

                // Correct velocity
                U -= rAU*fvc::grad(p);
                U.correctBoundaryConditions();
            }
        } while (++oCorr < nOuterCorr);


        // Correct density
        rho = rhoO + psi*p;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
