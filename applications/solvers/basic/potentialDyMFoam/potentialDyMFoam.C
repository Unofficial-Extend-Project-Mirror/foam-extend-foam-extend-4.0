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
    potentialDyMFoam

Description
    Transient solver for potential flow with dynamic mesh.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validOptions.insert("resetU", "");
    argList::validOptions.insert("writep", "");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "createFields.H"
#   include "initTotalVolume.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
#       include "readPISOControls.H"
#       include "checkTotalVolume.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;

        bool meshChanged = mesh.update();
        reduce(meshChanged, orOp<bool>());

        p.internalField() = 0;

        if (args.optionFound("resetU"))
        {
            U.internalField() = vector::zero;
        }

#       include "volContinuity.H"
#       include "meshCourantNo.H"

        // Solve potential flow equations
        adjustPhi(phi, U, p);

        for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
        {
            p.storePrevIter();

            fvScalarMatrix pEqn
            (
                fvm::laplacian
                (
                    dimensionedScalar
                    (
                        "1",
                        dimTime/p.dimensions()*dimensionSet(0, 2, -2, 0, 0),
                        1
                    ),
                    p
                )
             ==
                fvc::div(phi)
            );

            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve();

            if (nonOrth == nNonOrthCorr)
            {
                phi -= pEqn.flux();
            }
            else
            {
                p.relax();
            }
        }

        Info<< "continuity error = "
            << mag(fvc::div(phi))().weightedAverage(mesh.V()).value()
            << endl;

        U = fvc::reconstruct(phi);
        U.correctBoundaryConditions();

        Info<< "Interpolated U error = "
            << (sqrt(sum(sqr((fvc::interpolate(U) & mesh.Sf()) - phi)))
            /sum(mesh.magSf())).value()
            << endl;

        // Calculate velocity magnitude
        {
            volScalarField magU = mag(U);

            Info<< "mag(U): max: " << gMax(magU.internalField())
                << " min: " << gMin(magU.internalField()) << endl;
        }

        if (args.optionFound("writep"))
        {
            // Find reference patch
            label refPatch = -1;
            scalar maxMagU = 0;

            // Go through all velocity patches and find the one that fixes
            // velocity to the largest value

            forAll (U.boundaryField(), patchI)
            {
                const fvPatchVectorField& Upatch = U.boundaryField()[patchI];

                if (Upatch.fixesValue())
                {
                    // Calculate mean velocity
                    scalar u = sum(mag(Upatch));
                    label patchSize = Upatch.size();

                    reduce(u, sumOp<scalar>());
                    reduce(patchSize, sumOp<label>());

                    if (patchSize > 0)
                    {
                        scalar curMag = u/patchSize;

                        if (curMag > maxMagU)
                        {
                            refPatch = patchI;

                            maxMagU = curMag;
                        }
                    }
                }
            }

            if (refPatch > -1)
            {
                // Calculate reference pressure
                const fvPatchVectorField& Upatch = U.boundaryField()[refPatch];
                const fvPatchScalarField& pPatch = p.boundaryField()[refPatch];

                scalar patchE = sum(mag(pPatch + 0.5*magSqr(Upatch)));
                label patchSize = Upatch.size();

                reduce(patchE, sumOp<scalar>());
                reduce(patchSize, sumOp<label>());

                scalar e = patchE/patchSize;

                Info<< "Using reference patch " << refPatch
                    << " with mag(U) = " << maxMagU
                    << " p + 0.5*U^2 = " << e << endl;

                p.internalField() = e - 0.5*magSqr(U.internalField());
                p.correctBoundaryConditions();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
