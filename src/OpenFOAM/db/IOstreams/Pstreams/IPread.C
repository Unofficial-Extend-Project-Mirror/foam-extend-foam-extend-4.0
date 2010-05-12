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
    Read token and binary block from IPstream

\*---------------------------------------------------------------------------*/

#include "IPstream.H"
#include "error.H"
#include "int.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private member functions  * * * * * * * * * * * //

inline void IPstream::checkEof()
{
    if (bufPosition_ == messageSize_)
    {
        setEof();
    }
}


template<class T>
inline void IPstream::readFromBuffer(T& t)
{
    const size_t align = sizeof(T);
    bufPosition_ = align + ((bufPosition_ - 1) & ~(align - 1));

    t = reinterpret_cast<T&>(buf_[bufPosition_]);
    bufPosition_ += sizeof(T);
    checkEof();
}


inline void IPstream::readFromBuffer(void* data, size_t count, size_t align)
{
    if (align > 1)
    {
        bufPosition_ = align + ((bufPosition_ - 1) & ~(align - 1));
    }

    register const char* bufPtr = &buf_[bufPosition_];

    register char* dataPtr = reinterpret_cast<char*>(data);
    register size_t i = count;
    while (i--) *dataPtr++ = *bufPtr++;
    bufPosition_ += count;
    checkEof();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

IPstream::~IPstream()
{}


Istream& IPstream::read(char& c)
{
    c = buf_[bufPosition_];
    bufPosition_++;
    checkEof();
    return *this;
}


Istream& IPstream::read(word& w)
{
    size_t ws;
    readFromBuffer(ws);
    w = &buf_[bufPosition_];
    bufPosition_ += ws + 1;
    checkEof();
    return *this;
}


Istream& IPstream::read(string& s)
{
    size_t ss;
    readFromBuffer(ss);
    s = &buf_[bufPosition_];
    bufPosition_ += ss + 1;
    checkEof();
    return *this;
}


Istream& IPstream::read(label& l)
{
    readFromBuffer(l);
    return *this;
}


Istream& IPstream::read(floatScalar& s)
{
    readFromBuffer(s);
    return *this;
}


Istream& IPstream::read(doubleScalar& s)
{
    readFromBuffer(s);
    return *this;
}


Istream& IPstream::read(char* data, std::streamsize count)
{
    if (format() != BINARY)
    {
        FatalErrorIn("IPstream::read(char*, std::streamsize)")
            << "stream format not binary"
            << Foam::abort(FatalError);
    }

    readFromBuffer(data, count, 8);
    return *this;
}


Istream& IPstream::rewind()
{
    bufPosition_ = 0;
    return *this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
