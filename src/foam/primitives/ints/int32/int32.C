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

#include "int32.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const int32_t Foam::pTraits<int32_t>::zero = 0;
const int32_t Foam::pTraits<int32_t>::one = 1;
const int32_t Foam::pTraits<int32_t>::min = INT32_MIN;
const int32_t Foam::pTraits<int32_t>::max = INT32_MAX;
const int32_t Foam::pTraits<int32_t>::rootMin = pTraits<int32_t>::min;
const int32_t Foam::pTraits<int32_t>::rootMax = pTraits<int32_t>::max;

const char* Foam::pTraits<int32_t>::componentNames[] = { "x" };

Foam::pTraits<int32_t>::pTraits(const int32_t& p)
:
    p_(p)
{}

Foam::pTraits<int32_t>::pTraits(Istream& is)
{
    is >> p_;
}


// ************************************************************************* //
