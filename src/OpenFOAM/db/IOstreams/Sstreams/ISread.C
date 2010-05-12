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
    Read token and binary block from ISstream

\*---------------------------------------------------------------------------*/

#include "ISstream.H"
#include "int.H"
#include "token.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Istream& ISstream::read(char& c)
{
    c = nextValid();
    return *this;
}


Istream& ISstream::read(word& w)
{
    static const int MAX_WORD = 1024;
    static char wordBuffer[MAX_WORD];

    register int i = 0;
    register int bc = 0;
    char c;

    while (get(c) && word::valid(c))
    {
        if (fail())
        {
            if (i < MAX_WORD-1)
            {
                wordBuffer[i] = '\0';
            }
            else
            {
                wordBuffer[MAX_WORD-1] = '\0';
            }

            FatalIOErrorIn("ISstream::read(word&)", *this)
                << "problem while reading word '" << wordBuffer << '\''
                << exit(FatalIOError);

            return *this;
        }

        if (i >= MAX_WORD)
        {
            wordBuffer[MAX_WORD-1] = '\0';

            FatalIOErrorIn("ISstream::read(word&)", *this)
                << "word " << wordBuffer << " too long" << endl
                << "    maximum number of characters allowed = " << MAX_WORD
                << exit(FatalIOError);

            return *this;
        }

        if (c == token::BEGIN_LIST)
        {
            bc++;
        }
        else if (c == token::END_LIST)
        {
            bc--;

            if (bc == -1)
            {
                break;
            }
        }

        wordBuffer[i++] = c;
    }

    if (i == 0)
    {
        FatalIOErrorIn("ISstream::read(word&)", *this)
            << "invalid first character found : " << c
            << exit(FatalIOError);
    }

    wordBuffer[i] = '\0';        // Terminator.
    w = wordBuffer;
    putback(c);

    return *this;
}


Istream& ISstream::read(string& s)
{
    static const int MAX_STR = 1024;
    static const int MAX_ERROR_STR = 80;
    static char stringBuffer[MAX_STR];

    char c;

    if (!get(c))
    {
        stringBuffer[0] = '\0';

        FatalIOErrorIn("ISstream::read(string& s)", *this)
            << "cannot read start of string"
            << exit(FatalIOError);

        return *this;
    }

    if (c != token::BEGIN_STRING)
    {
        stringBuffer[0] = '\0';

        FatalIOErrorIn("ISstream::read(string& s)", *this)
            << "Incorrect start of string character"
            << exit(FatalIOError);

        return *this;
    }


    register int i = 0;
    bool escape = false;

    while (get(c))
    {
        if (c == token::END_STRING && !escape)
        {
            stringBuffer[i] = '\0';
            s = stringBuffer;
            return *this;
        }

        if (c == '\n' && !escape)
        {
            stringBuffer[i] = '\0';
            stringBuffer[MAX_ERROR_STR] = '\0';

            FatalIOErrorIn("ISstream::read(string& s)", *this)
                << "found a '\\n' while reading string \""
                << stringBuffer << '"'
                << exit(FatalIOError);

            return *this;
        }

        if (c != '\\')
        {
            escape = false;

            stringBuffer[i] = c;

            if (i++ == MAX_STR)
            {
                stringBuffer[MAX_STR - 1] = '\0';

                FatalIOErrorIn("ISstream::read(string& s)", *this)
                    << " string " << stringBuffer << "is too long" << endl
                    << "    maximum number of characters allowed = " << MAX_STR
                    << exit(FatalIOError);

                return *this;
            }
        }
        else
        {
            escape = true;
        }
    }

    stringBuffer[i] = '\0';
    stringBuffer[MAX_ERROR_STR] = '\0';

    FatalIOErrorIn("ISstream::read(string& s)", *this)
        << "problem while reading string \"" << stringBuffer << "...\""
        << exit(FatalIOError);

    return *this;
}


Istream& ISstream::read(label& l)
{
    is_ >> l;
    setState(is_.rdstate());
    return *this;
}


Istream& ISstream::read(floatScalar& s)
{
    is_ >> s;
    setState(is_.rdstate());
    return *this;
}


Istream& ISstream::read(doubleScalar& s)
{
    is_ >> s;
    setState(is_.rdstate());
    return *this;
}


// read binary block
Istream& ISstream::read(char* buf, std::streamsize count)
{
    if (format() != BINARY)
    {
        FatalIOErrorIn("ISstream::read(char*, std::streamsize)", *this)
            << "stream format not binary"
            << exit(FatalIOError);
    }

    readBegin("binaryBlock");
    is_.read(buf, count);
    readEnd("binaryBlock");

    setState(is_.rdstate());

    return *this;
}


//- Rewind the ISstream so that it may be read again
Istream& ISstream::rewind()
{
    stream().rdbuf()->pubseekpos(0);

    return *this;
}


// Set flags of output stream
ios_base::fmtflags ISstream::flags() const
{
    return is_.flags();
}


// Set flags of given field of output stream
ios_base::fmtflags ISstream::flags(const ios_base::fmtflags f)
{
    return is_.flags(f);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
