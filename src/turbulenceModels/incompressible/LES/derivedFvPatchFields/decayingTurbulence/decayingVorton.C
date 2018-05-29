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

\*---------------------------------------------------------------------------*/

#include "decayingVorton.H"
#include "Random.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::Random Foam::decayingVorton::ranGen(0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decayingVorton::decayingVorton
(
    const scalar length,
    const vector location,
    const vector velocity,
    const scalar xmax
)
:
    length_(length),
    location_(location),
    velocity_(velocity),
    xmax_(xmax),
    omega_(vector::zero)
{
    omega_ = 2*ranGen.vector01() - vector::one;
    omega_ /= mag(omega_);
}


Foam::decayingVorton::decayingVorton(Istream& s)
:
    length_(readScalar(s)),
    location_(s),
    velocity_(s),
    xmax_(readScalar(s)),
    omega_(s)
{}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

bool Foam::operator==(const decayingVorton& a, const decayingVorton& b)
{
    return
        mag(a.length_ - b.length_) < SMALL
     && mag(a.location_ - b.location_) < SMALL
     && mag(a.velocity_ - b.velocity_) < SMALL
     && mag(a.xmax_ - b.xmax_) < SMALL
     && mag(a.omega_ - b.omega_) < SMALL;
}


bool Foam::operator!=(const decayingVorton& a, const decayingVorton& b)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const decayingVorton& vt)
{
    os  << vt.length_ << nl
        << vt.location_ << nl
        << vt.velocity_ << nl
        << vt.xmax_ << nl
        << vt.omega_ << endl;

    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, decayingVorton& vt)
{
    is  >> vt.length_
        >> vt.location_
        >> vt.velocity_
        >> vt.xmax_
        >> vt.omega_;

    return is;
}


// ************************************************************************* //

