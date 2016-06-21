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

#include "thermoParcelInjectionData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::thermoParcelInjectionData::thermoParcelInjectionData(Istream& is)
:
    kinematicParcelInjectionData(is)
{
    is.check("reading T");
    is >> T_;

    is.check("reading cp");
    is >> cp_;

    is.check("thermoParcelInjectionData(Istream& is)");
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const thermoParcelInjectionData& data
)
{
    os << static_cast<const kinematicParcelInjectionData&>(data);

    os << data.T_ << data.cp_;

    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, thermoParcelInjectionData& data)
{
    is >> static_cast<kinematicParcelInjectionData&>(data);

    is.check("reading T");
    is >> data.T_;

    is.check("reading cp");
    is >> data.cp_;

    is.check("operator(Istream&, thermoParcelInjectionData&)");

    return is;
}


// ************************************************************************* //
