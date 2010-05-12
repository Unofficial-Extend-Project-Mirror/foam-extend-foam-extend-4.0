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

#include "ISstream.H"
#include "token.H"
#include <cctype>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Istream& ISstream::read(token& t)
{
    static char numberBuffer[100];

    // Return the put back token if it exists
    if (Istream::getBack(t))
    {
        return *this;
    }

    // Assume that the streams supplied are in working order.
    // Lines are counted by '\n'

    // Get next 'valid character': i.e. proceed through any white space
    // and/or comments until a semantically valid character is hit upon.

    char c = nextValid();

    // Set the line number of this token to the current stream line number
    t.lineNumber() = lineNumber();

    // return on error
    if (!c)
    {
        t.setBad();
        return *this;
    }

    // Analyse input starting with this character.
    switch (c)
    {
        // First check for punctuation characters.

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
     // case token::SUBTRACT : // Handled later as the posible start of a number
        case token::MULTIPLY :
        case token::DIVIDE :
        {
            t = token::punctuationToken(c);
            return *this;
        }

        // Strings: enclosed by double quotes.
        case token::BEGIN_STRING :
        {
            putback(c);
            string* sPtr = new string;

            if (!read(*sPtr).bad())
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

        // Numbers: do not distinguish at this point between Types.
        case '-' :
        case '.' :
        case '0' : case '1' : case '2' : case '3' : case '4' :
        case '5' : case '6' : case '7' : case '8' : case '9' :
        {
            bool isScalar = false;

            if (c == '.')
            {
                isScalar = true;
            }

            int i=0;
            numberBuffer[i++] = c;

            while
            (
                is_.get(c)
             && (
                    isdigit(c)
                 || c == '.'
                 || c == 'e' 
                 || c == 'E'
                 || c == '+'
                 || c == '-'
                )
            )
            {
                numberBuffer[i++] = c;

                if (!isdigit(c))
                {
                    isScalar = true;
                }
            }
            numberBuffer[i] = '\0';

            setState(is_.rdstate());

            if (!is_.bad())
            {
                is_.putback(c);

                if (i == 1 && numberBuffer[0] == '-')
                {
                    t = token::punctuationToken(token::SUBTRACT);
                }
                else if (isScalar)
                {
                    t = scalar(atof(numberBuffer));
                }
                else
                {
                    t = label(atol(numberBuffer));
                }
            }
            else
            {
                t.setBad();
            }

            return *this;
        }

        // Should be a word (which can be a single character)
        default:
        {
            putback(c);
            word* wPtr = new word;

            if (!read(*wPtr).bad())
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
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
