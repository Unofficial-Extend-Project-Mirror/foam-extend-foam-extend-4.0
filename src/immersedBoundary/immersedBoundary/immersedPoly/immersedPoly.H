/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    immersedPoly

Description
    immersedPoly holds basic data for immersed face/cell objects
    Wet side is "inside", indicated by negative distance to surface

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
            WET = 0,           // Inside
            DRY = 1            // Outside
        };

        //- Number of iterations in the iterative intersection
        static const debug::optimisationSwitch nIter_;

        //- Tolerance for considering a point on surface. Used as a fraction
        //  of edge length
        static const debug::tolerancesSwitch tolerance_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
