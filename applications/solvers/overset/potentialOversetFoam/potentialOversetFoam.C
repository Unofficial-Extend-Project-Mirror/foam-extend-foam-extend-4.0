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
    potentialOversetFoam

Description
    Potential flow solver with overset mesh support.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"

#include "oversetMesh.H"
#include "oversetFvPatchFields.H"
#include "oversetAdjustPhi.H"
#include "globalOversetAdjustPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validOptions.insert("resetU", "");
    argList::validOptions.insert("writep", "");

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    simpleControl simple(mesh);

#   include "createOversetMasks.H"
#   include "createFields.H"

#   include "writeOversetMasks.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Calculating potential flow" << endl;

    // Do correctors over the complete set
    while (simple.correctNonOrthogonal())
    {
        phi = (linearInterpolate(U) & mesh.Sf());

        // Adjust overset fluxes
        oversetAdjustPhi(phi, U); // Fringe flux adjustment
        globalOversetAdjustPhi(phi, U, p); // Global flux adjustment

        Info<< "Initial flux contour continuity error = "
            << mag(sum(phi.boundaryField()))
            << endl;

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

        // Adjust non-orthogonal fringe fluxes if necessary
        om.correctNonOrthoFluxes(pEqn, U);

        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve();

        // Correct the flux
        phi -= pEqn.flux();

        volScalarField divFlux
        (
            "divFlux",
            cellOversetMask*fvc::div(phi)
        );
        divFlux.write();

        // Perform overset interpolation (after flux reconstruction)
        oversetFvPatchScalarField::oversetInterpolate(p);

        if (!simple.finalNonOrthogonalIter())
        {
            p.relax();
        }

        Info<< "p min " << gMin(p.internalField())
            << " max " << gMax(p.internalField())
            << " masked min "
            << gMin(cellOversetMask.internalField()*p.internalField())
            << " max "
            << gMax(cellOversetMask.internalField()*p.internalField())
            << endl;

        Info<< "continuity error = "
            << mag
               (
                    cellOversetMask*fvc::div(phi)
               )().weightedAverage(mesh.V()).value()
            << endl;

        Info<< "Contour continuity error = "
            << mag(sum(phi.boundaryField()))
            << endl;

        U = fvc::reconstruct(phi);
        U.correctBoundaryConditions();
        oversetFvPatchVectorField::oversetInterpolate(U); // Overset update
        // Note: if implicit fringe is switched on, we are doing the
        // interpolation twice (once in correctBoundaryConditions and once
        // in oversetInterpolate). Reorganize. VV, 4/Oct/2016.

        Info<< "Interpolated U error = "
            << (
                   sqrt
                   (
                       sum
                       (
                           sqr
                           (
                               faceOversetMask*
                               (
                                   (fvc::interpolate(U) & mesh.Sf())
                                 - phi
                               )
                           )
                       )
                   )/sum(mesh.magSf())
               ).value()
            << endl;
    }

    // Calculate velocity magnitude
    {
        volScalarField magU = cellOversetMask*mag(U);

        Info << "Masked mag(U): max: " << gMax(magU.internalField())
            << " min: " << gMin(magU.internalField()) << endl;
    }

    // Force the write
    U.write();
    p.write();
    phi.write();
    cellOversetMask.write();

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
            oversetFvPatchScalarField::oversetInterpolate(p); // Overset update
            // Note: if implicit fringe is switched on, we are doing the
            // interpolation twice (once in correctBoundaryConditions and once
            // in oversetInterpolate). Reorganize. VV, 4/Oct/2016.
        }
        else
        {
            Info<< "No reference patch found.  Writing potential function"
                << endl;
        }

        p.write();
    }

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
