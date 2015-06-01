/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "IPstream.H"
#include "token.H"
#include <cctype>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Istream& IPstream::read(token& t)
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
            word* wPtr = new word;
            if (read(*wPtr))
            {
                if (token::compound::isCompound(*wPtr))
                {
                    t = token::compound::New(*wPtr, *this).ptr();
                    delete wPtr;
                }
                else
                {
                    t = wPtr;
                }
            }
            else
            {
                delete wPtr;
                t.setBad();
            }
            return *this;
        }

        // String
        case token::STRING :
        {
            string* sPtr = new string;
            if (read(*sPtr))
            {
                t = sPtr;
            }
            else
            {
                delete sPtr;
                t.setBad();
            }
            return *this;
        }

        // Label
        case token::LABEL :
        {
            label l;
            if (read(l))
            {
                t = l;
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
            floatScalar s;
            if (read(s))
            {
                t = s;
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
            doubleScalar s;
            if (read(s))
            {
                t = s;
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
