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

#include "quaternion.H"
#include "IOstreams.H"
#include "OStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::quaternion::typeName = "quaternion";
const Foam::quaternion Foam::quaternion::zero(0, vector::zero);
const Foam::quaternion Foam::quaternion::I(1, vector::zero);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::quaternion::quaternion(Istream& is)
{
    operator>>(is, *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::name(const quaternion& q)
{
    OStringStream buf;
    buf << '(' << q.w() << ',' << q.v() << ')';
    return buf.str();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, quaternion& q)
{
    // Read beginning of quaternion
    is.readBegin("quaternion");

    is  >> q.w() >> q.v();

    // Read end of quaternion
    is.readEnd("quaternion");

    // Check state of Istream
    is.check("operator>>(Istream&, quaternion&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const quaternion& q)
{
    os  << token::BEGIN_LIST
        << q.w() << token::SPACE << q.v()
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
