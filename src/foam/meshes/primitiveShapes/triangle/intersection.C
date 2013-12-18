/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "intersection.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::intersection::planarTol_
(
    debug::tolerances("intersectionPlanarTol", 0.2)
);


const Foam::scalar Foam::intersection::missTol_
(
    debug::tolerances("intersectionMissTol", SMALL)
);


template<>
const char* Foam::NamedEnum<Foam::intersection::direction, 2>::names[] =
{
    "vector",
    "contactSphere"
};


const Foam::NamedEnum<Foam::intersection::direction, 2>
Foam::intersection::directionNames_;


template<>
const char* Foam::NamedEnum<Foam::intersection::algorithm, 3>::names[] =
{
    "fullRay",
    "halfRay",
    "visible"
};


const Foam::NamedEnum<Foam::intersection::algorithm, 3>
Foam::intersection::algorithmNames_;



// ************************************************************************* //
