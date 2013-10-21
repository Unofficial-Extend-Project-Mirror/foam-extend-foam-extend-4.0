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
    splineInterpolateXY

Description
    1-D cubic spline interpolation of scalar x/y-data
    For details of the algorithm, see "Numerical Recipes in C" chapter 3.3

SourceFiles
    linearWaveVelocityForwardSpeedFvPatchField.C

Author
    Dominik Christ, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef splineInterpolateXY_H
#define splineInterpolateXY_H

#include "scalarField.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class splineInterpolateXY Declaration
\*---------------------------------------------------------------------------*/

class splineInterpolateXY
{
public:

    //- Types of spline boundary condition
    enum splineBCType
    {
        NOT_A_KNOT,
        NATURAL
    };

    static const NamedEnum<splineBCType, 2> splineBCTypeNames;


private:

    // Private data

        //- Given x-data
        scalarField X_;
        
        //- Given y-data
        scalarField Y_;

        //- Second derivative data
        scalarField Y2_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        splineInterpolateXY(const splineInterpolateXY&);

        //- Disallow default bitwise assignment
        void operator=(const splineInterpolateXY&);


        //- Initialize spline with given values
        void setData
        (
            const splineBCType& startBC,
            const splineBCType& endBC
        );


public:

    // Constructors

        //- Construct from components
        splineInterpolateXY
        (
            const scalarField& x,
            const scalarField& y,
            const word& startSplineBCs,
            const word& endSplineBCs
        );

    // Destructor - default


    // Member Functions

        //- Return interpolated value at given location xStar
        scalar interpolate(const scalar xStar) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
