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

#include "reactingMultiphaseParcelInjectionData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::reactingMultiphaseParcelInjectionData::
reactingMultiphaseParcelInjectionData(Istream& is)
:
    reactingParcelInjectionData(is)
{
    is.check("reading YGas's");
    is >> YGas_;

    is.check("reading YLiquid's");
    is >> YLiquid_;

    is.check("reading YSolid's");
    is >> YSolid_;

    is.check("reactingMultiphaseParcelInjectionData(Istream& is)");
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const reactingMultiphaseParcelInjectionData& data
)
{
    os << static_cast<const reactingParcelInjectionData&>(data);

    os << data.YGas_ << data.YLiquid_ << data.YSolid_;

    return os;
}


Foam::Istream& Foam::operator>>
(
    Istream& is,
    reactingMultiphaseParcelInjectionData& data
)
{
    is >> static_cast<reactingParcelInjectionData&>(data);

    is.check("reading YGas's");
    is >> data.YGas_;

    is.check("reading YLiquid's");
    is >> data.YLiquid_;

    is.check("reading YSolid's");
    is >> data.YSolid_;

    is.check("operator(Istream&, reactingMultiphaseParcelInjectionData&)");

    return is;
}


// ************************************************************************* //
