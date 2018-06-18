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
    RichardsFoam

Description
    Transient solver for flow in unsaturated porous media
    With chord slope formulation of the Richards equation.
    van Genuchten laws for unsaturated hydraulic properties parametrisation
    Global computation of the convergence criterium
    Adaptative time stepping with a stabilisation procedure
    NB 1: use backward scheme for time discretisation
    NB 2: use only mesh with constant cell volumes

References
   version 0.0 (develloped with OpenFOAM 2.0.1)
   Details may be found in:
   Orgogozo, L., Renon, N., Soulaine, C., Hénon, F., Tomer, S.K., Labat, D.,
    Pokrovsky, O.S., Sekhar, M., Ababou, R., Quintard, M., Submitted.
   Mechanistic modelling of water fluxes at the watershed scale: An open source
    massively parallel solver for Richards equation.
   Submitted to Computer Physics Communications.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

//    pimpleControl pimple(mesh);

#   include "readPicardControls.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"
#   include "createTimeControls.H"

    Info<< "\nStarting time loop\n" << endl;

    // starting of the time loop.
    while (runTime.loop())
    {
        // time step control operations.
#       include "readTimeControls.H"
#       include "setDeltaT.H"

//        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Beginning of the stabilisation loop for the stabilised adaptive time
        // step procedure.
        for (int cyc = 0; cyc < nMaxCycle; cyc++)
        {
            // Beginning of the Picard loop.
            for (int pic = 0; pic < nIterPicard; pic++)
            {
#               include "psiEqn.H"
            }

            // Exit test for the loop associated with the stabilisation cycles
            // for the adaptive time step procedure.
            if (crit < precPicard)
            {
                break;
            }
            else
            {
                Info << "Criterion not reached, restart time loop iteration"
                    << "with a smaller time step / Error = " << crit
                    << nl << endl;

                runTime.setDeltaT((1/tFact)*runTime.deltaTValue());

                Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
            }
            // End of the stabilisation cycles loop.
        }

        // Warning test in case of convergence failure of the Picard loop.
        if (crit >= precPicard)
        {
            Info<< "Convergence failure / Error = " << crit << nl << endl;
            currentPicard = nIterPicard;
        }

        // Final updating of the result fields before going to the next time
        // iteration.
        psi_tmp = psi;

        thtil_tmp = 0.5*
        (
            (1 + sign(psi_tmp)) + (1 - sign(psi_tmp))*
            pow((1 + pow(mag(alpha*psi_tmp),n)), - (1 - (1/n)))
        );

        theta = (thetas - thetar)*thtil + thetar;

        U = - Krel*((fvc::grad(psi)) + vuz);

        // Writting  of the result.
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        // end of the time loop.
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
