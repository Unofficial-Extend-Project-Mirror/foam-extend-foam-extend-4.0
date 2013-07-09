/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "viewFactorRadiation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viewFactorRadiation::viewFactorRadiation
(
    const word& name,
    const dictionary& dict,
    const label size
)
:
    name_(name),
    Tinf_(readScalar(dict.lookup("Tinf"))),
    F_("F", dict, size),
    epsilon_(readScalar(dict.lookup("epsilon")))
{}


void Foam::viewFactorRadiation::write(Ostream& os) const
{
    os  << indent << name_ << nl
        << indent << token::BEGIN_BLOCK << incrIndent << nl;

    os.writeKeyword("Tinf") << Tinf_ << token::END_STATEMENT << nl;
    F_.writeEntry("F", os);
    os.writeKeyword("epsilon") << epsilon_ << token::END_STATEMENT << nl;

    os   << decrIndent << indent << token::END_BLOCK << endl;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const viewFactorRadiation& r)
{
    r.write(os);
    os.check("Ostream& operator<<(Ostream&, const viewFactorRadiation&)");
    return os;
}


// ************************************************************************* //
