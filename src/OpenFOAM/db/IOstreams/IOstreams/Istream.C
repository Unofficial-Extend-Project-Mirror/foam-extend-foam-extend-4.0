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

#include "Istream.H"
#include "bool.H"
#include "token.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set t to the put back token if there is one and return true, 
// otherwise return false
bool Istream::getBack(token& t)
{
    if (bad())
    {
        FatalIOErrorIn("void Istream::getBack(token& t)", *this)
            << "Attempt to get back from bad stream"
            << exit(FatalIOError);

        return false;
    }
    else if (putBack_)
    {
        t = putBackToken_;
        putBack_ = false;
        return true;
    }
    else
    {
        return false;
    }
}


// Keep the put back token
void Istream::putBack(const token& t)
{
    if (bad())
    {
        FatalIOErrorIn("void Istream::putBack(const token& t)", *this)
            << "Attempt to put back onto bad stream"
            << exit(FatalIOError);
    }
    else if (putBack_)
    {
        FatalIOErrorIn("void Istream::putBack(const token& t)", *this)
            << "Attempt to put back another token"
            << exit(FatalIOError);
    }
    else
    {
        putBackToken_ = t;
        putBack_ = true;
    }
}


// Functions for reading object delimiters ( ... )

Istream& Istream::readBegin(const char* funcName)
{
    token delimiter(*this);
    if (delimiter != token::BEGIN_LIST)
    {
        setBad();
        FatalIOErrorIn("Istream::readBegin(const char*)", *this)
            << "Expected a " << '\'' << token:: BEGIN_LIST << '\''
            << " while reading " << funcName
            << ", found " << delimiter.info()
            << exit(FatalIOError);
    }

    return *this;
}


Istream& Istream::readEnd(const char* funcName)
{
    token delimiter(*this);
    if (delimiter != token::END_LIST)
    {
        setBad();
        FatalIOErrorIn("Istream::readEnd(const char*)", *this)
            << "Expected a " << '\'' << token::END_LIST << '\''
            << " while reading " << funcName
            << ", found " << delimiter.info()
            << exit(FatalIOError);
    }

    return *this;
}


Istream& Istream::readEndBegin(const char* funcName)
{
    readEnd(funcName);
    return readBegin(funcName);
}


// Functions for reading List delimiters ( ... ) or { ... }

char Istream::readBeginList(const char* funcName)
{
    token delimiter(*this);

    if (delimiter != token::BEGIN_LIST && delimiter != token::BEGIN_BLOCK)
    {
        setBad();
        FatalIOErrorIn("Istream::readBeginList(const char*)", *this) 
            << "Expected a " << '\'' << token::BEGIN_LIST << '\''
            << " or a " << '\'' << token::BEGIN_BLOCK << '\''
            << " while reading " << funcName
            << ", found " << delimiter.info()
            << exit(FatalIOError);

        return '\0';
    }

    return delimiter.pToken();
}


char Istream::readEndList(const char* funcName)
{
    token delimiter(*this);

    if (delimiter != token::END_LIST && delimiter != token::END_BLOCK)
    {
        setBad();
        FatalIOErrorIn("Istream::readEndList(const char*)", *this)
            << "Expected a " << '\'' << token::END_LIST << '\''
            << " or a " << '\'' << token::END_BLOCK << '\''
            << " while reading " << funcName
            << ", found " << delimiter.info()
            << exit(FatalIOError);

        return '\0';
    }

    return delimiter.pToken();
}


//- Return a non-const reference to const Istream
//  Needed for read-constructors where the stream argument is temporary:
//  e.g. thing thisThing(IFstream("thingFileName")());
Istream& Istream::operator()() const
{
    if (!good())
    {
        check("Istream::operator()");
        FatalIOError.exit();
    }

    return const_cast<Istream&>(*this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
