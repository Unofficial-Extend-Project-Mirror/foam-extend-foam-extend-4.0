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

Description
    Return location of contact sphere on the face

\*---------------------------------------------------------------------------*/

#include "face.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::face::contactSphereDiameter
(
    const point& p,
    const vector& n,
    const pointField& meshPoints
) const
{
    vector n1 = n/(Foam::mag(n) + SMALL);
    vector n2 = normal(meshPoints);

    n2 /= Foam::mag(n2) + SMALL;

    return 2*((centre(meshPoints) - p) & n2)/((n1 & n2) - 1.0);
}


// ************************************************************************* //
