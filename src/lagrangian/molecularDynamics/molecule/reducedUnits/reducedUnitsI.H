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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::reducedUnits::refLength() const
{
    return refLength_;
}


inline Foam::scalar Foam::reducedUnits::refTime() const
{
    return refTime_;
}


inline Foam::scalar Foam::reducedUnits::refMass() const
{
    return refMass_;
}


inline Foam::scalar Foam::reducedUnits::refTemp() const
{
    return refTemp_;
}


inline Foam::scalar Foam::reducedUnits::refEnergy() const
{
    return refEnergy_;
}


inline Foam::scalar Foam::reducedUnits::refForce() const
{
    return refForce_;
}


inline Foam::scalar Foam::reducedUnits::refVelocity() const
{
    return refVelocity_;
}


inline Foam::scalar Foam::reducedUnits::refVolume() const
{
    return refVolume_;
}


inline Foam::scalar Foam::reducedUnits::refPressure() const
{
    return refPressure_;
}


inline Foam::scalar Foam::reducedUnits::refMassDensity() const
{
    return refMassDensity_;
}


inline Foam::scalar Foam::reducedUnits::refNumberDensity() const
{
    return refNumberDensity_;
}


// ************************************************************************* //
