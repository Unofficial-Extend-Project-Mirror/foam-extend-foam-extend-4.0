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
    Write primitive and binary block from OPstream

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "OPstream.H"
#include "int.H"
#include "token.H"

#include <cctype>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

OPstream::OPstream
(
    const commsTypes commsType,
    const int toProcNo,
    const label bufSize,
    streamFormat format,
    versionNumber version
)
:
    Pstream(commsType, bufSize),
    Ostream(format, version),
    toProcNo_(toProcNo)
{
    setOpened();
    setGood();

    if (!bufSize)
    {
        buf_.setSize(1000);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Ostream& OPstream::write(const token&)
{
    notImplemented("Ostream& OPstream::write(const token&)");
    setBad();
    return *this;
}


Ostream& OPstream::write(const char c)
{
    if (!isspace(c))
    {
        writeToBuffer(c);
    }

    return *this;
}


Ostream& OPstream::write(const char* s)
{
    word nonWhiteChars(string::validate<word>(s));

    if (nonWhiteChars.size() == 0)
    {
        return *this;
    }
    else if (nonWhiteChars.size() == 1)
    {
        return write(nonWhiteChars.c_str()[1]);
    }
    else
    {
        return write(nonWhiteChars);
    }
}


Ostream& OPstream::write(const word& w)
{
    write(char(token::WORD));

    size_t ws = w.size();
    writeToBuffer(ws);
    writeToBuffer(w.c_str(), ws + 1, 1);

    return *this;
}


Ostream& OPstream::write(const string& s)
{
    write(char(token::STRING));

    size_t ss = s.size();
    writeToBuffer(ss);
    writeToBuffer(s.c_str(), ss + 1, 1);

    return *this;
}


Ostream& OPstream::write(const label l)
{
    write(char(token::LABEL));

    writeToBuffer(l);

    return *this;
}


Ostream& OPstream::write(const floatScalar s)
{
    write(char(token::FLOAT_SCALAR));

    writeToBuffer(s);

    return *this;
}


Ostream& OPstream::write(const doubleScalar s)
{
    write(char(token::DOUBLE_SCALAR));

    writeToBuffer(s);

    return *this;
}


Ostream& OPstream::write(const char* data, std::streamsize count)
{
    if (format() != BINARY)
    {
        FatalErrorIn("Ostream::write(const char*, std::streamsize)")
            << "stream format not binary"
            << Foam::abort(FatalError);
    }

    writeToBuffer(data, count, 8);

    return *this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
