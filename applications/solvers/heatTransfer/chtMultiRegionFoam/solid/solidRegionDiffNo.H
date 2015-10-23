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
    Calculates and outputs the mean and maximum Diffusion Numbers for the solid
    regions

\*---------------------------------------------------------------------------*/

#ifndef solidRegionDiff_H
#define solidRegionDiff_H

#include "fvMesh.H"

namespace Foam
{
    scalar solidRegionDiffNo
    (
        const fvMesh& mesh,
        const Time& runTime,
        const volScalarField& Cprho,
        const volScalarField& Kappa
    );
}

#endif

// ************************************************************************* //
