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

#include "ListOps.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::invert
(
    const label len,
    const UList<label>& map
)
{
    labelList inverse(len, -1);

    forAll(map, i)
    {
        label newPos = map[i];

        if (newPos >= 0)
        {
            if (inverse[newPos] >= 0)
            {
                FatalErrorIn("invert(const label, const UList<label>&)")
                    << "Map is not one-to-one. At index " << i
                    << " element " << newPos << " has already occurred before"
                    << nl << "Please use invertOneToMany instead"
                    << abort(FatalError);
            }

            inverse[newPos] = i;
        }
    }
    return inverse;
}


Foam::labelListList Foam::invertOneToMany
(
    const label len,
    const UList<label>& map
)
{
    labelList nElems(len, 0);

    forAll(map, i)
    {
        if (map[i] >= 0)
        {
            nElems[map[i]]++;
        }
    }

    labelListList inverse(len);

    forAll(nElems, i)
    {
        inverse[i].setSize(nElems[i]);
        nElems[i] = 0;
    }

    forAll(map, i)
    {
        label newI = map[i];

        if (newI >= 0)
        {
            inverse[newI][nElems[newI]++] = i;
        }
    }

    return inverse;
}


Foam::labelList Foam::identity(const label len)
{
    labelList map(len);

    forAll(map, i)
    {
        map[i] = i;
    }
    return map;
}


// ************************************************************************* //
