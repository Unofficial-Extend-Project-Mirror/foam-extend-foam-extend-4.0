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

\*----------------------------------------------------------------------------*/

#include "referredMolecule.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::referredMolecule::referredMolecule()
{}


Foam::referredMolecule::referredMolecule
(
    const label id,
    const vector& position,
    const List<vector>& sitePositions
)
:
    id_(id),
    position_(position),
    sitePositions_(sitePositions)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::referredMolecule::~referredMolecule()
{}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is,
    referredMolecule& rM
)
{
    is >> rM.id_ >> rM.position_ >> rM.sitePositions_;

    is.check("Istream& operator<<(Istream& f, const referredMolecule& sRL");

    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const referredMolecule& rM
)
{
    os  << rM.id()
        << token::SPACE << rM.position()
        << token::SPACE << rM.sitePositions();

    os.check("Ostream& operator<<(Ostream& f, const referredMolecule& rM");

    return os;
}


// ************************************************************************* //
