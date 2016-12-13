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

\*---------------------------------------------------------------------------*/

#include "dimensionSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::dimensionSet Foam::dimless(0, 0, 0, 0, 0, 0, 0);

const Foam::dimensionSet Foam::dimMass(1, 0, 0, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimLength(0, 1, 0, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimTime(0, 0, 1, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimTemperature(0, 0, 0, 1, 0, 0, 0);
const Foam::dimensionSet Foam::dimMoles(0, 0, 0, 0, 1, 0, 0);
const Foam::dimensionSet Foam::dimCurrent(0, 0, 0, 0, 0, 1, 0);
const Foam::dimensionSet Foam::dimLuminousIntensity(0, 0, 0, 0, 0, 0, 1);

const Foam::dimensionSet Foam::dimArea(0, 2, 0, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimVolume(0, 3, 0, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimVol(0, 3, 0, 0, 0, 0, 0);

const Foam::dimensionSet Foam::dimDensity(1, -3, 0, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimForce(1, 1, -2, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimEnergy(1, 2, -2, 0, 0, 0, 0);

const Foam::dimensionSet Foam::dimVelocity(0, 1, -1, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimAcceleration(0, 1, -2, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimPressure(1, -1, -2, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimGasConstant(0, 2, -2, -1, 0, 0, 0);
const Foam::dimensionSet Foam::dimSpecificHeatCapacity(0, 2, -2, -1, 0, 0, 0);
const Foam::dimensionSet Foam::dimThermalConductivity(1, 1, -3, -1, 0, 0, 0);

const Foam::dimensionSet Foam::dimVoltage(1, 2, -3, 0, 0, -1, 0);


// ************************************************************************* //
