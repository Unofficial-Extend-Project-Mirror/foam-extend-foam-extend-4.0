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

Description
    A simple wrapper around bool so that it can
    be read as on/off, yes/no or y/n.

\*---------------------------------------------------------------------------*/

#include "Switch.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Switch::Switch(Istream& is)
:
    wordSwitch_(is),
    logicalSwitch_(boolValue(wordSwitch_))
{}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Istream& operator>>(Istream& is, Switch& s)
{
    is >> s.wordSwitch_;
    s.logicalSwitch_ = s.boolValue(s.wordSwitch_);

    is.check("Istream& operator>>(Istream& is, Switch& s)");

    return is;
}


Ostream& operator<<(Ostream& os, const Switch& s)
{
    os << s.wordSwitch_;

    os.check("Ostream& operator<<(Ostream& os, const Switch& s)");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
