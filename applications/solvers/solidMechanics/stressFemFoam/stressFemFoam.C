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
    stressFemFoam

Description
    Steady-state segregated finite-element solver of linear-elastic,
    small-strain deformation of a solid body.

    Simple linear elasticity structural analysis code.
    Solves for the displacement vector field U, also generating the
    stress tensor field sigma.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "pointMesh.H"
#include "pointFields.H"
#include "tetPolyMesh.H"
#include "tetPointFields.H"
#include "tetFem.H"
#include "tetFec.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

#   include "createFemMesh.H"
#   include "createFemFields.H"

#   include "readMechanicalProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nCalculating displacement field\n" << endl;

    while(runTime.loop())
    {
        Info<< "Iteration: " << runTime.timeName() << nl << endl;

        tetFemVectorMatrix UEqn
        (
            tetFem::laplacian(mu, U)
          + tetFem::laplacianTranspose(mu, U)
          + tetFem::laplacianTrace(lambda, U)
        );

        UEqn.solve();

#       include "calculateStress.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
