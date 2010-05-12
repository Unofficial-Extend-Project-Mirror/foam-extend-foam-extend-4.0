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
    Reads an bool from an input stream, for a given version
    number and File format. If an ascii File is being read,
    then the line numbers are counted and an erroneous read
    ised.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "bool.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Istream& operator>>(Istream& is, bool& b)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isLabel())
    {
        b = bool(t.labelToken());
    }
    else if (t.isWord())
    {
        if (t.wordToken() == "true" || t.wordToken() == "on")
        {
            b = true;
        }
        else if (t.wordToken() == "false" || t.wordToken() == "off")
        {
            b = false;
        }
        else
        {
            is.setBad();
            FatalIOErrorIn("operator>>(Istream&, bool&)", is)
                << "expected 'true' or 'false', found " << t.wordToken()
                << exit(FatalIOError);

            return is;
        }
    }
    else
    {
        is.setBad();
        FatalIOErrorIn("operator>>(Istream&, bool&)", is)
            << "wrong token type - expected bool found " << t
            << exit(FatalIOError);

        return is;
    }



    // Check state of Istream
    is.check("Istream& operator>>(Istream&, bool&)");

    return is;
}


Ostream& operator<<(Ostream& os, const bool b)
{
    os.write(label(b));
    os.check("Ostream& operator<<(Ostream&, const bool&)");
    return os;
}


bool readBool(Istream& is)
{
    bool rb;
    is >> rb;

    return rb;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
