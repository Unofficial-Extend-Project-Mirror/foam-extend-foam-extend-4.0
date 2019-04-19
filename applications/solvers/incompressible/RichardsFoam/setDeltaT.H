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

Global
    setDeltaT

Description
    Reset the timestep to maintain a low number of Picard iteration and a good
    convergence.  Reduction of time-step is immediate, but increase is damped to
    avoid unstable oscillations.

\*---------------------------------------------------------------------------*/

scalar deltaTFact = 1.;

if (adjustTimeStep)
{
    deltaTFact = 1.;

    // If There is a too high number of Picard iteration, decrease of the time
    // step.
    if (currentPicard > nIterPicard-2)
    {
        deltaTFact = 1./tFact;
    }

    if (currentPicard >= 3)
    {
        sc = 0;
    }

    if ((currentPicard < 3) && (sc < stabilisationThreshold))
    {
        sc = sc+1;
    }

    // If there is more than 'stabilisationThreshold' iteration that the number
    // of Picard iteration is low, increase of the time step.
    if ((currentPicard < 3) && (sc == stabilisationThreshold))
    {
        deltaTFact = tFact;
        sc = 0;
    }

    // Reset of the time step if needed (time step is always lower than
    // 'maxDeltaT').
    if (deltaTFact != 1.)
    {
        runTime.setDeltaT
        (
            min
            (
                deltaTFact*runTime.deltaTValue(),
                maxDeltaT
            )
        );
    }

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}

// ************************************************************************* //
