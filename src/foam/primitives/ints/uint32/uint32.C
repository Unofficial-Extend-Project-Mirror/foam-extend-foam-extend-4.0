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

#include "uint32.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const uint32_t Foam::pTraits<uint32_t>::zero = 0;
const uint32_t Foam::pTraits<uint32_t>::one = 1;
const uint32_t Foam::pTraits<uint32_t>::min = INT32_MIN;
const uint32_t Foam::pTraits<uint32_t>::max = INT32_MAX;
const uint32_t Foam::pTraits<uint32_t>::rootMin = pTraits<uint32_t>::min;
const uint32_t Foam::pTraits<uint32_t>::rootMax = pTraits<uint32_t>::max;

const char* Foam::pTraits<uint32_t>::componentNames[] = { "x" };

Foam::pTraits<uint32_t>::pTraits(const uint32_t& p)
:
    p_(p)
{}

Foam::pTraits<uint32_t>::pTraits(Istream& is)
{
    is >> p_;
}


// ************************************************************************* //
