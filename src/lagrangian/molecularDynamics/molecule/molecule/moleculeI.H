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

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline molecule::molecule
(
    const Cloud<molecule>& c,
    const vector& position,
    const label celli,
    const scalar mass,
    const vector& U,
    const vector& A,
    const vector& tetherPosition,
    const label tethered,
    const label id
)
:
    Particle<molecule>(c, position, celli),
    mass_(mass),
    U_(U),
    A_(A),
    tetherPosition_(tetherPosition),
    potentialEnergy_(0.0),
    rf_(tensor::zero),
    tethered_(tethered),
    id_(id)
{}


inline molecule::trackData::trackData
(
    moleculeCloud& molCloud,
    label part
)
:
    Particle<molecule>::trackData(refCast<Cloud<molecule> >(molCloud)),
    molCloud_(molCloud),
    part_(part)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline label molecule::id() const
{
    return id_;
}


inline scalar molecule::mass() const
{
    return mass_;
}


inline const vector& molecule::U() const
{
    return U_;
}


inline vector& molecule::U()
{
    return U_;
}


inline const vector& molecule::A() const
{
    return A_;
}


inline vector& molecule::A()
{
    return A_;
}


inline scalar molecule::potentialEnergy() const
{
    return potentialEnergy_;
}


inline scalar& molecule::potentialEnergy()
{
    return potentialEnergy_;
}


inline const tensor& molecule::rf() const
{
    return rf_;
}


inline tensor& molecule::rf()
{
    return rf_;
}


inline label molecule::tethered() const
{
    return tethered_;
}


inline const vector& molecule::tetherPosition() const
{
    return tetherPosition_;
}


inline vector& molecule::tetherPosition()
{
    return tetherPosition_;
}


inline moleculeCloud& molecule::trackData::molCloud()
{
    return molCloud_;
}


inline label molecule::trackData::part() const
{
    return part_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
