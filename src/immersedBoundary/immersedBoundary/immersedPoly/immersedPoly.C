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

#include "immersedPoly.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::debug::optimisationSwitch
Foam::immersedPoly::nIter_
(
    "immersedPolyNIter",
    3
);


const Foam::debug::tolerancesSwitch
Foam::immersedPoly::tolerance_
(
    "immersedPolyTolerance",
    1e-4
);


const Foam::debug::tolerancesSwitch
Foam::immersedPoly::liveFactor_
(
    "immersedPolyLiveFactor",
    0.01
);


const Foam::debug::tolerancesSwitch
Foam::immersedPoly::badCutFactor_
(
    "immersedPolyBadCutFactor",
    0.01
);


// ************************************************************************* //
