/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Global
    temperatureAndPressure.H

Description
    Accumulates values for temperature and pressure measurement, and
    calculates and outputs the average values at output times.
    Requires temperatureAndPressureVariables.H to be declared before the
    timeloop.

\*---------------------------------------------------------------------------*/

accumulatedTotalMomentum += singleStepTotalMomentum;

accumulatedTotalMass += singleStepTotalMass;

accumulatedTotalKE += singleStepTotalKE;

accumulatedTotalPE += singleStepTotalPE;

accumulatedTotalrDotfSum += singleStepTotalrDotf;

accumulatedNMols += singleStepNMols;

if (runTime.outputTime())
{
    // calculate averages

    if (accumulatedNMols)
    {
        Info << "calculating averages" << endl;

        averageTemperature =
        (
            2.0/(3.0 * moleculeCloud::kb * accumulatedNMols)
            *
            (
                accumulatedTotalKE
                -
                0.5*magSqr(accumulatedTotalMomentum)/accumulatedTotalMass
            )
        );

        averagePressure =
        (
            (
                (accumulatedNMols/nAveragingSteps)
                *
                moleculeCloud::kb * averageTemperature
                +
                accumulatedTotalrDotfSum/(6.0 * nAveragingSteps)
            )
            /
            meshVolume
        );

        // output values

        Info << "----------------------------------------" << nl
            << "Averaged properties" << nl
            << "Average |velocity| = "
                << mag(accumulatedTotalMomentum)/accumulatedTotalMass
                << " m/s" << nl
            << "Average temperature = "
                << averageTemperature << " K" << nl
            << "accumulatedTotalrDotfSum = "
                << accumulatedTotalrDotfSum << nl
            << "Average pressure = "
                << averagePressure << " N/m^2" << nl
            << "----------------------------------------" << endl;
    }
    else
    {
        Info << "Not averaging temperature and pressure: "
            << "no molecules in system" << endl;
    }

    // reset counters

    accumulatedTotalMomentum = vector::zero;

    accumulatedTotalMass = 0.0;

    accumulatedTotalKE = 0.0;

    accumulatedTotalPE = 0.0;

    accumulatedTotalrDotfSum = 0.0;

    accumulatedNMols = 0;
}

// ************************************************************************* //
