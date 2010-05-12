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
    Reads an ulong from an input stream, for a given version
    number and File format. If an ascii File is being read,
    then the line numbers are counted and an erroneous read
    ised.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "ulong.H"
#include "IOstreams.H"

#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

word name(const unsigned long i)
{
    std::ostringstream osBuffer;
    osBuffer << i;
    return osBuffer.str();
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Istream& operator>>(Istream& is, unsigned long& i)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isLabel())
    {
        i = static_cast<unsigned long>(t.labelToken());
    }
    else
    {
        is.setBad();
        FatalIOErrorIn("operator>>(Istream&, unsigned long&)", is)
            << "wrong token type - expected unsigned long found " << t.info()
            << exit(FatalIOError);

        return is;
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, unsigned long&)");

    return is;
}


unsigned long readUlong(Istream& is)
{
    unsigned long ri;
    is >> ri;

    return ri;
}


Ostream& operator<<(Ostream& os, const unsigned long i)
{
    os.write(label(i));
    os.check("Ostream& operator<<(Ostream&, const unsigned long)");
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
