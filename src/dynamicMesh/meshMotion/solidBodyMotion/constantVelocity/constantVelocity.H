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
    Foam::solidBodyMotionFunction::constantVelocity

Description
    Prescribed constant translational and angular velocity.

SourceFiles
    constantVelocity.C

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef constantVelocity_H
#define constantVelocity_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                         Class constantVelocity Declaration
\*---------------------------------------------------------------------------*/

class constantVelocity
:
    public solidBodyMotionFunction
{
    // Private data

        //- Centre of gravity
        point origin_;

        //- Translational velocity vector
        vector transVelocity_;

        //- Rotational velocity vector
        vector rotVelocity_;

        //- Start motion time
        scalar startMotionTime_;

        //- Is the rotational velocity given in degrees/sec
        Switch inDegrees_;


    // Private Member Functions

        //- Disallow copy construct
        constantVelocity(const constantVelocity&);

        //- Disallow default bitwise assignment
        void operator=(const constantVelocity&);

        //- Calculate tranformation
        septernion calcTransformation(const scalar t) const;


public:

    //- Runtime type information
    TypeName("constantVelocity");


    // Constructors

        //- Construct from components
        constantVelocity
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );


    //- Destructor
    virtual ~constantVelocity()
    {}


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Return the solid-body motion velocity
        virtual septernion velocity() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
