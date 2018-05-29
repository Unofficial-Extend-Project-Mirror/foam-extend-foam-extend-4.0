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
    sonicDyMFoam

Description
    Transient solver for trans-sonic/supersonic for laminar or turbulent
    flow of a compressible gas with support for mesh motion and
    topological changes.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient simulations.  The pressure-energy coupling is done
    using the Rusche manoeuvre (isentropic compression/expansion).

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

    Updated from sonicFoamAutoMotion by Hrvoje Jasak

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "specie.H"
#include "basicPsiThermo.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createDynamicFvMesh.H"

    pimpleControl pimple(mesh);

#   include "createFields.H"
#   include "initContinuityErrs.H"
#   include "createControls.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readControls.H"
#       include "readFieldBounds.H"
#       include "compressibleCourantNo.H"
#       include "setDeltaT.H"

        runTime++;
        Info<< "deltaT = " << runTime.deltaT().value() << nl << endl;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        bool meshChanged = mesh.update();

#       include "volContinuity.H"

        if (checkMeshCourantNo)
        {
#           include "meshCourantNo.H"
        }

        // Mesh motion update
        if (meshChanged)
        {
            T.max(TMin);
            p.max(pMin);
            e == max(e, thermo.Cv()*TMin);

            thermo.correct();
            rho = thermo.rho();

            if (correctPhi)
            {
// #           include "correctPhi.H"
            }
        }

        if (meshChanged)
        {
#           include "compressibleCourantNo.H"
        }

        // --- PIMPLE loop
        while (pimple.loop())
        {
#           include "rhoEqn.H"
#           include "eEqn.H"
#           include "UEqn.H"

            // --- PISO loop
            volScalarField rUA = 1.0/UEqn.A();

            surfaceScalarField psisf = fvc::interpolate(psis);
            surfaceScalarField rhof = fvc::interpolate(rho);

            // Needs to be outside of loop since p is changing,
            // but psi and rho are not
            surfaceScalarField rhoReff = rhof - psisf*fvc::interpolate(p);

            while (pimple.correct())
            {
#               include "pEqn.H"
            }

            turbulence->correct();
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
