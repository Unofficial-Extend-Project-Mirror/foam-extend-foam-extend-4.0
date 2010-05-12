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
    axisCoordinateRotation

Description
    The 3 rotations are defined in the axis convention
    (around X, around Y and around Z).
    For reference and illustration, see
    http://mathworld.wolfram.com/RotationMatrix.html

AUTHOR
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef axisCoordinateRotation_H
#define axisCoordinateRotation_H

#include "coordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class axisCoordinateRotation Declaration
\*---------------------------------------------------------------------------*/

class axisCoordinateRotation
:
    public coordinateRotation
{
    // Private member functions

        //- Calculate transformation tensor
        void calcTransform
        (
            const scalar phiAngle,
            const scalar thetaAngle,
            const scalar psiAngle,
            const bool inDegrees = true
        );


public:

    //- Runtime type information
    TypeName("axisRotation");


    // Constructors

        //- Construct null
        axisCoordinateRotation();

        //- Construct from components of rotation vector
        axisCoordinateRotation
        (
            const scalar phiAngle,
            const scalar thetaAngle,
            const scalar psiAngle,
            const bool inDegrees = true
        );

        //- Construct from dictionary
        axisCoordinateRotation(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
