/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "IPstream.H"
#include "int.H"
#include "token.H"
#include <cctype>


// * * * * * * * * * * * * * Private member functions  * * * * * * * * * * * //

inline void Foam::IPstream::checkEof()
{
    if (bufPosition_ == messageSize_)
    {
        setEof();
    }
}


template<class T>
inline void Foam::IPstream::readFromBuffer(T& t)
{
    const size_t align = sizeof(T);
    bufPosition_ = align + ((bufPosition_ - 1) & ~(align - 1));

    t = reinterpret_cast<T&>(buf_[bufPosition_]);
    bufPosition_ += sizeof(T);
    checkEof();
}


inline void Foam::IPstream::readFromBuffer
(
    void* data,
    size_t count,
    size_t align
)
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


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::IPstream::~IPstream()
{}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Istream& Foam::IPstream::read(token& t)
{
    // Return the put back token if it exists
    if (Istream::getBack(t))
    {
        return *this;
    }

    char c;

    // return on error
    if (!read(c))
    {
        t.setBad();
        return *this;
    }

    // Set the line number of this token to the current stream line number
    t.lineNumber() = lineNumber();

    // Analyse input starting with this character.
    switch (c)
    {
        // Punctuation
        case token::END_STATEMENT :
        case token::BEGIN_LIST :
        case token::END_LIST :
        case token::BEGIN_SQR :
        case token::END_SQR :
        case token::BEGIN_BLOCK :
        case token::END_BLOCK :
        case token::COLON :
        case token::COMMA :
        case token::ASSIGN :
        case token::ADD :
        case token::SUBTRACT :
        case token::MULTIPLY :
        case token::DIVIDE :
        {
            t = token::punctuationToken(c);
            return *this;
        }

        // Word
        case token::WORD :
        {
            word* pval = new word;
            if (read(*pval))
            {
                if (token::compound::isCompound(*pval))
                {
                    t = token::compound::New(*pval, *this).ptr();
                    delete pval;
                }
                else
                {
                    t = pval;
                }
            }
            else
            {
                delete pval;
                t.setBad();
            }
            return *this;
        }

        // String
        case token::STRING :
        {
            string* pval = new string;
            if (read(*pval))
            {
                t = pval;
            }
            else
            {
                delete pval;
                t.setBad();
            }
            return *this;
        }

        // Label
        case token::LABEL :
        {
            label val;
            if (read(val))
            {
                t = val;
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // floatScalar
        case token::FLOAT_SCALAR :
        {
            floatScalar val;
            if (read(val))
            {
                t = val;
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // doubleScalar
        case token::DOUBLE_SCALAR :
        {
            doubleScalar val;
            if (read(val))
            {
                t = val;
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // longDoubleScalar
        case token::LONG_DOUBLE_SCALAR :
        {
            longDoubleScalar val;
            if (read(val))
            {
                t = val;
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // Character (returned as a single character word) or error
        default:
        {
            if (isalpha(c))
            {
                t = word(c);
                return *this;
            }

            setBad();
            t.setBad();

            return *this;
        }
    }
}


Foam::Istream& Foam::IPstream::read(char& c)
{
    c = buf_[bufPosition_];
    bufPosition_++;
    checkEof();
    return *this;
}


Foam::Istream& Foam::IPstream::read(word& str)
{
    size_t len;
    readFromBuffer(len);
    str = &buf_[bufPosition_];
    bufPosition_ += len + 1;
    checkEof();
    return *this;
}


Foam::Istream& Foam::IPstream::read(string& str)
{
    size_t len;
    readFromBuffer(len);
    str = &buf_[bufPosition_];
    bufPosition_ += len + 1;
    checkEof();
    return *this;
}


Foam::Istream& Foam::IPstream::read(label& val)
{
    readFromBuffer(val);
    return *this;
}


Foam::Istream& Foam::IPstream::read(floatScalar& val)
{
    readFromBuffer(val);
    return *this;
}


Foam::Istream& Foam::IPstream::read(doubleScalar& val)
{
    readFromBuffer(val);
    return *this;
}


Foam::Istream& Foam::IPstream::read(longDoubleScalar& val)
{
    readFromBuffer(val);
    return *this;
}


Foam::Istream& Foam::IPstream::read(char* data, std::streamsize count)
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


Foam::Istream& Foam::IPstream::rewind()
{
    bufPosition_ = 0;
    return *this;
}


void Foam::IPstream::print(Ostream& os) const
{
    os  << "Reading from processor " << fromProcNo_
        << " using communicator " << comm_
        <<  " and tag " << tag_
        << Foam::endl;
}


// ************************************************************************* //
