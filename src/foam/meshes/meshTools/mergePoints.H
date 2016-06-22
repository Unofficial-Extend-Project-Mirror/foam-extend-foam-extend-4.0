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
    Merge points. See below.

SourceFiles
    mergePoints.C

\*---------------------------------------------------------------------------*/

#ifndef mergePoints_H
#define mergePoints_H

#include "scalarField.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Function mergePoints Declaration
\*---------------------------------------------------------------------------*/

//- Sort & merge points. All points closer than/equal mergeTol get merged.
//  Outputs the new unique points and a map from old to new. Returns
//  true if anything merged, false otherwise.
bool mergePoints
(
    const UList<point>& points,
    const scalar mergeTol,
    const bool verbose,
    labelList& pointMap,
    List<point>& newPoints,
    const point& origin = point::zero
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
