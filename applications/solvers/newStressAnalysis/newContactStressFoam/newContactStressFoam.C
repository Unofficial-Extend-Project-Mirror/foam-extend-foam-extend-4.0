/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    stressedFoam

Description
    Transient/steady-state solver for solid bodies in contact.

    Solves for the displacement vector field U, also generating the
    stress tensor field sigma.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rheologyModel.H"
#include "contactProblem.H"
#include "componentReferenceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating displacement field\n" << endl;

    volScalarField rho = rheology.rho();

    // Force n-sqaured projection
//     polyPatch::setNSquaredProjection(true);

    while (runTime.loop())
    {
        Info<< "Iteration: " << runTime.timeName() << nl << endl;

#       include "readStressedFoamControls.H"

        volScalarField mu = rheology.mu();
        volScalarField lambda = rheology.lambda();

        int iCorr=0;
        scalar initialResidual=0;

        contact.correct();

        do
        {
            fvVectorMatrix UEqn
            (
                fvm::d2dt2(rho, U)
              ==
                fvm::laplacian(2*mu + lambda, U, "laplacian(DU,U)")

              + fvc::div
                (
                    mu*gradU.T() + lambda*(I*tr(gradU)) - (mu + lambda)*gradU,
                    "div(sigma)"
                )
            );

#           include "setComponentReference.H"

            initialResidual = UEqn.solve().initialResidual();

            gradU = fvc::grad(U);

#           include "calculateSigma.H"

            rheology.correct();

            rho = rheology.rho();
            mu = rheology.mu();
            lambda = rheology.lambda();

        } while (initialResidual > convergenceTolerance && ++iCorr < nCorr);

#       include "calculateStress.H"
#       include "calculateContactArea.H"

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
