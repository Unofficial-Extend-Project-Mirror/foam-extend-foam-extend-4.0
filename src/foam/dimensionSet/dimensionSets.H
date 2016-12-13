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

Global
    dimensionSets

Description
    Useful dimension sets

SourceFiles
    dimensionSets.C

\*---------------------------------------------------------------------------*/


#ifndef dimensionSets_H
#define dimensionSets_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

extern const dimensionSet dimless;

extern const dimensionSet dimMass;
extern const dimensionSet dimLength;
extern const dimensionSet dimTime;
extern const dimensionSet dimTemperature;
extern const dimensionSet dimMoles;
extern const dimensionSet dimCurrent;
extern const dimensionSet dimLuminousIntensity;

extern const dimensionSet dimArea;
extern const dimensionSet dimVolume;
extern const dimensionSet dimVol;

extern const dimensionSet dimDensity;
extern const dimensionSet dimForce;
extern const dimensionSet dimEnergy;

extern const dimensionSet dimVelocity;
extern const dimensionSet dimAcceleration;
extern const dimensionSet dimPressure;
extern const dimensionSet dimGasConstant;
extern const dimensionSet dimSpecificHeatCapacity;
extern const dimensionSet dimThermalConductivity;

extern const dimensionSet dimVoltage;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
