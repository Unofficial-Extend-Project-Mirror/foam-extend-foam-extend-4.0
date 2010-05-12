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

#include "error.H"
#include "OSstream.H"
#include "token.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Ostream& OSstream::write(const token&)
{
    return *this;
}


Ostream& OSstream::write(const char c)
{
    os_ << c;
    if (c == token::NL)
    {
        lineNumber_++;
    }
    setState(os_.rdstate());
    return *this;
}


Ostream& OSstream::write(const char* s)
{
    lineNumber_ += string(s).count(token::NL);
    os_ << s;
    setState(os_.rdstate());
    return *this;
}


Ostream& OSstream::write(const word& w)
{
    os_ << w;
    setState(os_.rdstate());
    return *this;
}


Ostream& OSstream::write(const string& s)
{
    os_ << '\"';

    for
    (
        string::const_iterator iter = s.begin();
        iter != s.end();
        ++iter
    )
    {
        register char c = *iter;

        if (c == token::NL)
        {
            os_ << '\\';
            lineNumber_++;
        }

        if (c == '"')
        {
            os_ << '\\';
        }

        if (c != '\\')
        {
            os_ << c;
        }
    }

    os_ << '\"';

    setState(os_.rdstate());
    return *this;
}


Ostream& OSstream::write(const label l)
{
    os_ << l;
    setState(os_.rdstate());
    return *this;
}


Ostream& OSstream::write(const floatScalar s)
{
    os_ << s;
    setState(os_.rdstate());
    return *this;
}


Ostream& OSstream::write(const doubleScalar s)
{
    os_ << s;
    setState(os_.rdstate());
    return *this;
}


Ostream& OSstream::write(const char* buf, std::streamsize count)
{
    if (format() != BINARY)
    {
        FatalIOErrorIn("Ostream::write(const char*, std::streamsize)", *this)
            << "stream format not binary"
            << abort(FatalIOError);
    }

    os_ << token::BEGIN_LIST;
    os_.write(buf, count);
    os_ << token::END_LIST;

    setState(os_.rdstate());

    return *this;
}


//- Add indentation characters
void OSstream::indent()
{
    for (register unsigned short i = 0; i < indentLevel_*indentSize_; i++)
    {
        os_ << ' ';
    }
}


// Flush stream
void OSstream::flush()
{
    os_.flush();
}


// Add carriage return and flush stream
void OSstream::endl()
{
    write('\n');
    os_.flush();
}


// Set flags of output stream
ios_base::fmtflags OSstream::flags() const
{
    return os_.flags();
}


// Set flags of given field of output stream
ios_base::fmtflags OSstream::flags(const ios_base::fmtflags f)
{
    return os_.flags(f);
}


//- Get width of output field
int OSstream::width() const
{
    return os_.width();
}

// Set width of output field (and return old width)
int OSstream::width(const int w)
{
    return os_.width(w);
}

//- Get precision of output field
int OSstream::precision() const
{
    return os_.precision();
}

// Set precision of output field (and return old precision)
int OSstream::precision(const int p)
{
    return os_.precision(p);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
