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
    channelFoam

Description
    Incompressible LES solver for flow in a channel.
    Consistent formulation without time-step and relaxation dependence by Jasak

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "LESModel.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Random.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "readTransportProperties.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createGradP.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        sgsModel->correct();

        // Convection-diffusion matrix
        fvVectorMatrix HUEqn
        (
            fvm::div(phi, U)
          + sgsModel->divDevBeff()
         ==
            flowDirection*gradP
        );

        // Time derivative matrix
        fvVectorMatrix ddtUEqn(fvm::ddt(U));

        if (piso.momentumPredictor())
        {
            solve(ddtUEqn + HUEqn == -fvc::grad(p));
        }

        // Prepare clean Ap without time derivative contribution
        // HJ, 26/Oct/2015
        volScalarField aU = HUEqn.A();

        // --- PISO loop

        while (piso.correct())
        {
            U = HUEqn.H()/aU;
            phi = (fvc::interpolate(U) & mesh.Sf());

            adjustPhi(phi, U, p);

            while (piso.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(1/aU, p) == fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve
                (
                    mesh.solutionDict().solver(p.select(piso.finalInnerIter()))
                );

                if (piso.finalNonOrthogonalIter())
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

        // Correct driving force for a constant mass flow rate

        // Extract the velocity in the flow direction
        dimensionedScalar magUbarStar =
            (flowDirection & U)().weightedAverage(mesh.V());

        // Calculate the pressure gradient increment needed to
        // adjust the average flow-rate to the correct value
        dimensionedScalar gragPplus =
            (magUbar - magUbarStar)*aU.weightedAverage(mesh.V());

        U += gragPplus/aU*flowDirection;

        gradP += gragPplus;

        Info<< "Uncorrected Ubar = " << magUbarStar.value() << tab
            << "pressure gradient = " << gradP.value() << endl;

        runTime.write();

        #include "writeGradP.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
