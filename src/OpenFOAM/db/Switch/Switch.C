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
#include "error.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Switch Switch::lookupOrAddToDict
(
    const word& name,
    dictionary& dict,
    const Switch& defaultValue
)
{
    return dict.lookupOrAddDefault<Switch>(name, defaultValue);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

word Switch::wordValue(const bool l) const
{
    word w("off");

    if (l)
    {
        w = "on";
    }

    return w;
}


bool Switch::boolValue(const word& w) const
{
    bool l = true;

    if (w == "on" || w == "yes" || w == "y" || w == "true")
    {
        l = true;
    }
    else if (w == "off" || w == "no" || w == "n" || w == "false")
    {
        l = false;
    }
    else
    {
        FatalErrorIn("Switch::boolValue(const word& w) const")
            << "unknown switch word " << w
            << abort(FatalError);
    }

    return l;
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

bool Switch::readIfPresent(const word& name, const dictionary& dict)
{
    return dict.readIfPresent(name, logicalSwitch_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
