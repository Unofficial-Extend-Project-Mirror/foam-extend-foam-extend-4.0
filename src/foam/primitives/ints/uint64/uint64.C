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

#include "uint64.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const uint64_t Foam::pTraits<uint64_t>::zero = 0;
const uint64_t Foam::pTraits<uint64_t>::one = 1;
const uint64_t Foam::pTraits<uint64_t>::min = INT64_MIN;
const uint64_t Foam::pTraits<uint64_t>::max = INT64_MAX;
const uint64_t Foam::pTraits<uint64_t>::rootMin = pTraits<uint64_t>::min;
const uint64_t Foam::pTraits<uint64_t>::rootMax = pTraits<uint64_t>::max;

const char* Foam::pTraits<uint64_t>::componentNames[] = { "x" };

Foam::pTraits<uint64_t>::pTraits(const uint64_t& p)
:
    p_(p)
{}

Foam::pTraits<uint64_t>::pTraits(Istream& is)
{
    is >> p_;
}


// ************************************************************************* //
