/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if WM_LABEL_SIZE == 32
const char* const Foam::pTraits<int64_t>::typeName = "int64";
const char* const Foam::pTraits<int32_t>::typeName = "label";
#elif WM_LABEL_SIZE == 64
const char* const Foam::pTraits<int64_t>::typeName = "label";
const char* const Foam::pTraits<int32_t>::typeName = "int32";
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::pow(label a, label b)
{
    label ans = 1;
    for (label i=0; i<b; i++)
    {
        ans *= a;
    }

    #ifdef FULLDEBUG
    if (b < 0)
    {
        FatalErrorIn("pow(label a, label b)")
            << "negative value for b is not supported"
            << abort(FatalError);
    }
    #endif

    return ans;
}


Foam::label Foam::factorial(label n)
{
    static label factTable[13] =
    {
        1, 1, 2, 6, 24, 120, 720, 5040, 40320,
        362880, 3628800, 39916800, 479001600
    };

    #ifdef FULLDEBUG
    if (n > 12 && n < 0)
    {
        FatalErrorIn("factorial(label n)")
            << "n value out of range"
            << abort(FatalError);
    }
    #endif

    return factTable[n];
}


// ************************************************************************* //
