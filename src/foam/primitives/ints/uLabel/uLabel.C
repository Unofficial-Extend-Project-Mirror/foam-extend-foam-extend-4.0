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
#include "uLabel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if WM_LABEL_SIZE == 32
const char* const Foam::pTraits<uint64_t>::typeName = "uint64";
const char* const Foam::pTraits<uint32_t>::typeName = "uLabel";
#elif WM_LABEL_SIZE == 64
const char* const Foam::pTraits<uint64_t>::typeName = "uLabel";
const char* const Foam::pTraits<uint32_t>::typeName = "uint32";
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::uLabel Foam::pow(uLabel a, uLabel b)
{
    uLabel ans = 1;
    for (uLabel i=0; i<b; i++)
    {
        ans *= a;
    }

    return ans;
}


Foam::uLabel Foam::factorial(uLabel n)
{
    static uLabel factTable[13] =
    {
        1, 1, 2, 6, 24, 120, 720, 5040, 40320,
        362880, 3628800, 39916800, 479001600
    };

    #ifdef FULLDEBUG
    if (n > 12)
    {
        FatalErrorIn("factorial(uLabel n)")
            << "n value out of range"
            << abort(FatalError);
    }
    #endif

    return factTable[n];
}


// ************************************************************************* //
