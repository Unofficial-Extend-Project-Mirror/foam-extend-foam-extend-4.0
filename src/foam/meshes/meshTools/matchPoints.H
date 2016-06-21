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
    Determine correspondence between points. See below.

SourceFiles
    matchPoints.C

\*---------------------------------------------------------------------------*/

#ifndef matchPoints_H
#define matchPoints_H

#include "scalarField.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Function matchPoints Declaration
\*---------------------------------------------------------------------------*/

//- Determine correspondence between pointFields. Gets passed
//  local matching distance per point. Returns map from 0 to 1. Returns
//  true if found corresponding point in 1 for every point in 0; false
//  otherwise. Compares mag(pt - origin) and does proper comparison for points
//  with (almost) identical mag. pts0 can be subset of pts1.
bool matchPoints
(
    const UList<point>& pts0,
    const UList<point>& pts1,
    const UList<scalar>& matchDistance,
    const bool verbose,
    labelList& from0To1,
    const point& origin = point::zero
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
