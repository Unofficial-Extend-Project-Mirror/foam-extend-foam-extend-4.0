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
    Foam::intersection

Description
    Foam::intersection

SourceFiles
   intersection.C

\*---------------------------------------------------------------------------*/

#ifndef intersection_H
#define intersection_H

#include "scalar.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class intersection Declaration
\*---------------------------------------------------------------------------*/

class intersection
{
public:

    // Static data

        //- Relative tolerance for point in a plane.
        //  The tolerance is relative to the object size.
        //  Used to determine if a half-ray from a point in plane of triangle
        //  intersects the triangle or not.
        //  Adjusting the tolerance from the code does not allow it to be const
        //  HJ, 5/Jan/2008
        static const scalar planarTol_;

        //- Intersection miss tolerance, of the order of 1e-15
        static const scalar missTol_;


        enum direction
        {
            VECTOR,
            CONTACT_SPHERE
        };

        enum algorithm
        {
            FULL_RAY,        // Intersecting with full ray
            HALF_RAY,        // Intersecting with half ray
            VISIBLE          // Intersecting with the visible side
        };


    // Static member functions

        //- Direction names
        static const NamedEnum<direction, 2> directionNames_;

        //- Projection algorithm names
        static const NamedEnum<algorithm, 3> algorithmNames_;

        //- Set the planar tolerance, returning the previous value
        static scalar setPlanarTol(const scalar t)
        {
            if (t < -VSMALL)
            {
                FatalErrorIn
                (
                    "scalar intersection::setPlanarTol(const scalar t)"
                )   << "Negative planar tolerance.  This is not allowed."
                    << abort(FatalError);
            }

            // Ugly violence to allow tolerance to remain a const scalar
            // HJ, 7/Sep/2008
            scalar oldTol = planarTol_;
            scalar* ptPtr = const_cast<scalar*>(&planarTol_);
            (*ptPtr) = t;

            return oldTol;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
