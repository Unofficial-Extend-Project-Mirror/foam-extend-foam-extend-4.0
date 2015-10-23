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

Description

\*---------------------------------------------------------------------------*/

#include "physicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Avogadro constant (number of molecules in a mole)
const Foam::dimensionedScalar Foam::physicalConstant::A
(
    "A",
    dimensionSet(0, 0, 0, 0, -1, 0, 0),
    6.0221415e23
);


// Faraday's constant (charge of electron)
const Foam::dimensionedScalar Foam::physicalConstant::F
(
    "F",
    dimensionSet(0, 0, 1, 0, 0, 1, 0),
    9.6485309e4
);


// ************************************************************************* //
