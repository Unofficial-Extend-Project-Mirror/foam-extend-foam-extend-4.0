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
    consistentIcoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.
    Consistent formulation without time-step and relaxation dependence by Jasak

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

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

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readPISOControls.H"
#       include "CourantNo.H"

        // Convection-diffusion matrix
        fvVectorMatrix HUEqn
        (
            fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        // ddt matrix
        fvVectorMatrix ddtUEqn(fvm::ddt(U));

        solve(ddtUEqn + HUEqn == -fvc::grad(p));

        // Prepare clean Ap without time derivative contribution
        // HJ, 26/Oct/2015
        volScalarField aU = HUEqn.A();

        // --- PISO loop

        for (int corr = 0; corr < nCorr; corr++)
        {
            U = HUEqn.H()/aU;
            phi = (fvc::interpolate(U) & mesh.Sf());

            adjustPhi(phi, U, p);

            for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(1/aU, p) == fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (nonOrth == nNonOrthCorr)
                {
                    phi -= pEqn.flux();
                }
            }

#           include "continuityErrs.H"

            // Note: cannot call H(U) here because the velocity is not complete
            // HJ, 22/Jan/2016
            U = 1.0/(aU + ddtUEqn.A())*
                (
                    U*aU - fvc::grad(p) + ddtUEqn.H()
                );
            U.correctBoundaryConditions();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
