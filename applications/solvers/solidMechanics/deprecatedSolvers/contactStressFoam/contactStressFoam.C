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
    stressedFoam

Description
    Transient/steady-state solver of linear-elastic, small-strain deformation
    of solid bodies in contact.

    Simple linear elasticity structural analysis code.
    Solves for the displacement vector field U,
    also generating the stress tensor field sigma.

    Code can be compiled for simple stress analysis,
    and dynamic stress analysis.

    Additionally, special handling of very non-orthogonal meshes may be added.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fixedGradientFvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "contactPatchPair.H"

#include "pointFields.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating displacement field\n" << endl;

    volTensorField gradU = fvc::grad(U);

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Iteration: " << runTime.timeName() << nl << endl;

#       include "../stressedFoam/readStressedFoamControls.H"

        int iCorr=0;
        scalar initialResidual=0;

#       include "contactBoundaries.H"

        do
        {
#           include "tractionBoundaries.H"

            fvVectorMatrix UEqn
            (
#               ifdef Dynamic
                fvm::d2dt2(U)
              ==
#               endif

                fvm::laplacian(2*mu + lambda, U, "laplacian(DU,U)")

              + fvc::div
                (
                    mu*gradU.T() + lambda*(I*tr(gradU)) - (mu + lambda)*gradU,
                    "div(sigma)"
                )
            );

            UEqn.setComponentReference(6, 0, vector::Z, 0);

            initialResidual = UEqn.solve().initialResidual();

            gradU = fvc::grad(U);

        } while (initialResidual > convergenceTolerance && ++iCorr < nCorr);

#       include "calculateStress.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
