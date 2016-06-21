/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    fsiFoam

Description
    Finite volume fluid structure interaction solver based on partitioned
    approach and strong coupling.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "fluidSolidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "createSolidMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    fluidSolidInterface fsi(mesh, solidMesh);

    Info<< "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.value()
            << " (dt = " << runTime.deltaT().value() << ")" << nl << endl;

        fsi.initializeFields();

        fsi.updateInterpolator();

        scalar residualNorm = 0;

        if (fsi.predictor())
        {
            fsi.updateForce();

            fsi.stress().evolve();

            residualNorm =
                fsi.updateResidual();
        }

        do
        {
            fsi.outerCorr()++;

            fsi.updateDisplacement();

            fsi.moveFluidMesh();

            fsi.flow().evolve();

            fsi.updateForce();

            fsi.stress().evolve();

            residualNorm =
                fsi.updateResidual();
        }
        while
        (
            (residualNorm > fsi.outerCorrTolerance())
         && (fsi.outerCorr() < fsi.nOuterCorr())
        );

        fsi.flow().updateFields();
        fsi.stress().updateTotalFields();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
