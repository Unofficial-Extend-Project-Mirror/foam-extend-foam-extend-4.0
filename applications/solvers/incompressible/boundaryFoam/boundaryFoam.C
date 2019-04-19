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
    boundaryFoam

Description
    Steady-state solver for 1D turbulent flow, typically to generate boundary
    layer conditions at an inlet, for use in a simulation.

    Boundary layer code to calculate the U, k and epsilon distributions.
    Used to create inlet boundary conditions for experimental comparisons
    for which U and k have not been measured.
    Turbulence model is runtime selectable.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "wallFvPatch.H"
#include "makeGraph.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        fvVectorMatrix divR = turbulence->divDevReff();
        divR.source() = flowMask & divR.source();

        fvVectorMatrix UEqn
        (
            divR == gradP
        );

        UEqn.relax();

        UEqn.solve();


        // Correct driving force for a constant mass flow rate

        dimensionedVector UbarStar = flowMask & U.weightedAverage(mesh.V());

        U += (Ubar - UbarStar);
        gradP += (Ubar - UbarStar)/(1.0/UEqn.A())().weightedAverage(mesh.V());

        Info<< "Uncorrected Ubar = " << (flowDirection & UbarStar.value())<< tab
            << "pressure gradient = " << (flowDirection & gradP.value())
            << endl;

        turbulence->correct();

        if (runTime.outputTime())
        {
            volSymmTensorField R
            (
                IOobject
                (
                    "R",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                turbulence->R()
            );

            runTime.write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
