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

Class
    immersedPoly

Description
    immersedPoly holds basic data for immersed face/cell objects
    WET side is "inside", indicated by negative distance to surface.
    Positive distance to surface is DRY.

SourceFiles
    immersedPoly.C

Author
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#ifndef immersedPoly_H
#define immersedPoly_H

#include "optimisationSwitch.H"
#include "tolerancesSwitch.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class immersedPoly Declaration
\*---------------------------------------------------------------------------*/

class immersedPoly
{
public:

    // Static data

        //- Wet type enumeration
        enum immersedType
        {
            UNKNOWN = -1,      // Unknown
            WET = 0,           // Inside
            DRY = 1,           // Outside
            CUT = 2            // Cut
        };

        //- Number of iterations in the iterative intersection
        static const debug::optimisationSwitch nIter_;

        //- Tolerance for considering a point on surface.
        //  Used as a fraction of edge length
        static const debug::tolerancesSwitch tolerance_;

        //- Live factor
        //  Fraction of face/cell size within the IB for which the face/cell is
        //  considered to be live
        static const debug::tolerancesSwitch liveFactor_;

        //- Bad cut factor
        //  Fraction of face/cell size error that declares a bad cut
        static const debug::tolerancesSwitch badCutFactor_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
