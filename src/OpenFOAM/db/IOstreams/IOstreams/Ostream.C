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

\*---------------------------------------------------------------------------*/

#include "Ostream.H"
#include "token.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Decrememt the indent level
void Foam::Ostream::decrIndent()
{
    if (indentLevel_ == 0)
    {
        cerr<< "Ostream::decrIndent() : attempt to decrement 0 indent level"
            << std::endl;
    }
    else
    {
        indentLevel_--;
    }
}


// Write the keyword to the Ostream followed by appropriate indentation
Foam::Ostream& Foam::Ostream::writeKeyword(const Foam::word& keyword)
{
    indent();
    write(keyword);

    label nSpaces = max(entryIndentation_ - label(keyword.size()), 1);

    for (label i=0; i<nSpaces; i++)
    {
        write(char(token::SPACE));
    }

    return *this;
}


// ************************************************************************* //
