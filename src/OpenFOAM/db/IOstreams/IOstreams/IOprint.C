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
    Prints out a description of the IOstream to Serr.

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void IOstream::print(Ostream& os) const
{
    os  << "IOstream: " << "Version "  << version_ << ", format ";

    switch (format_)
    {
        case ASCII:
            os  << "ASCII";
        break;

        case BINARY:
            os  << "BINARY";
        break;
    }

    os  << ", line "       << lineNumber();

    if (opened())
    {
        os  << ", OPENED";
    }

    if (closed())
    {
        os  << ", CLOSED";
    }

    if (good())
    {
        os  << ", GOOD";
    }

    if (eof())
    {
        os  << ", EOF";
    }

    if (fail())
    {
        os  << ", FAIL";
    }

    if (bad())
    {
        os  << ", BAD";
    }

    os  << endl;
}


void IOstream::print(Ostream& os, const int streamState) const
{
    if (streamState == ios_base::goodbit)
    {
        os  << "ios_base::goodbit set : the last operation on stream succeeded"
            << endl;
    }
    else if (streamState & ios_base::badbit)
    {
        os  << "ios_base::badbit set : characters possibly lost"
            << endl;
    }
    else if (streamState & ios_base::failbit)
    {
        os  << "ios_base::failbit set : some kind of formatting error"
            << endl;
    }
    else if (streamState & ios_base::eofbit)
    {
        os  << "ios_base::eofbit set : at end of stream"
            << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
