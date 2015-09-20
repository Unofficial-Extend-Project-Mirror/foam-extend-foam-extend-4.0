/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "radiationConstants.H"
#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Stefan-Boltzmann constant (default in [J/(K4 m2 s)])
const Foam::debug::constantsSwitch
Foam::radiation::sigmaSB_
(
    "sigmaSB",
    5.670E-08,
    "Stefan-Boltzmann constant [J/(K4 m2 s)]"
);

const Foam::dimensionedScalar Foam::radiation::sigmaSB
(
    "sigmaSB",
    dimensionSet(1, 0, -3, -4, 0, 0, 0),
    sigmaSB_()
);

// ************************************************************************* //
