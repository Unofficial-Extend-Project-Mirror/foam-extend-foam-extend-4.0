/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "CaCO3.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CaCO3, 0);
    addToRunTimeSelectionTable(solid, CaCO3,);
    addToRunTimeSelectionTable(solid, CaCO3, Istream);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CaCO3::CaCO3()
:
    solid(2710, 850, 1.3, 0.0, 1.0)
{
    if (debug)
    {
        WarningIn("CaCO3::CaCO3()")
            << "Properties of CaCO3 need to be checked!!!"
            << endl;
    }
}


Foam::CaCO3::CaCO3(const solid& s)
:
    solid(s)
{}


Foam::CaCO3::CaCO3(Istream& is)
:
    solid(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CaCO3::writeData(Ostream& os) const
{
    solid::writeData(os);
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const CaCO3& s)
{
    s.writeData(os);
    return os;
}


// ************************************************************************* //
