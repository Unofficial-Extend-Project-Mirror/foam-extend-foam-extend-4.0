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
    Implementation of lexer. Return next semantically valid character in
    the given input stream, or zero if an error is encountered.

\*---------------------------------------------------------------------------*/

#include "ISstream.H"

#include <cctype>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char Foam::ISstream::nextValid()
{
    char c = 0;

    for (;;)
    {
        // Get next non-whitespace character
        while (get(c) && isspace(c))
        {}

        // Return if stream is bad
        if (bad() || isspace(c))
        {
            return 0;
        }

        // Is this the start of a C/C++ comment
        if (c == '/')
        {
            // If cannot get another charater, return this one
            if (!get(c))
            {
                return '/';
            }

            if (c == '/')  // This is the start of a C++ style one line comment
            {
                while (get(c) && c != '\n')
                {}
            }
            else if (c == '*')  // This is the start of a C style comment
            {
                for (;;)
                {
                    if (get(c) && c == '*')
                    {
                        if (get(c) && c == '/')
                        {
                            break;
                        }
                        else
                        {
                            putback(c);
                        }
                    }

                    if (!good())
                    {
                        return 0;
                    }
                }
            }
            else  // A lone '/' so return it.
            {
                putback(c);
                return '/';
            }
        }
        else  // c is a valid character so return it
        {
            return c;
        }
    }
}


// ************************************************************************* //
