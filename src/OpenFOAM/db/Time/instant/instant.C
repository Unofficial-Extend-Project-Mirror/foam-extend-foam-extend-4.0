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

#include "instantList.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::instant::typeName = "instant";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::instant::instant()
{}

Foam::instant::instant(const scalar tval, const word& tname)
:
    value_(tval),
    name_(tname)
{}

Foam::instant::instant(const scalar tval)
:
    value_(tval),
    name_(Time::timeName(tval))
{}

Foam::instant::instant(const word& tname)
:
    value_(atof(tname.c_str())),
    name_(tname)
{}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

int Foam::operator==(const instant& I1, const instant& I2)
{
    return
    (
        I1.value_ < I2.value_ + SMALL
     && I1.value_ > I2.value_ - SMALL
    );
}


int Foam::operator != (const instant& I1, const instant& I2)
{
    // Invert the '==' operator ('0'='false')
    return I1 == I2 ? 0 : 1;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, instant& I)
{
    is >> I.value_ >> I.name_;

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const instant& I)
{
   os << I.value_ << tab << I.name_;

   return os;
}


// ************************************************************************* //
