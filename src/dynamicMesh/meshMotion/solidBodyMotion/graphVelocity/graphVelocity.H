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

Class
    Foam::solidBodyMotionFunction::graphVelocity

Description
    Prescribed translational and rotational velocity given graphs for surge,
    sway, heave, roll, pitch and yaw. Hence, the velocity is given in graphs,
    not the motion (trajectory).

SourceFiles
    graphVelocity.C

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef graphVelocity_H
#define graphVelocity_H

#include "solidBodyMotionFunction.H"
#include "graph.H"
#include "primitiveFields.H"
#include "point.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                         Class graphVelocity Declaration
\*---------------------------------------------------------------------------*/

class graphVelocity
:
    public solidBodyMotionFunction
{
    // Private data

        //- Centre of gravity
        point origin_;

        //- Integrated translation
        mutable vector translation_;

        //- Integrated rotation
        mutable vector rotation_;

        //- Integrated old translation
        mutable vector translationOld_;

        //- Integrated old rotation
        mutable vector rotationOld_;

        //- Time index, avoid integrating more than once per time step
        mutable label localTimeIndex_;


        // Prescribed translation members

            // - Graph for surge translational velocity (x-translation)
            graph surge_;

            // - Graph for sway translational velocity (y-translation)
            graph sway_;

            // - Graph for heave translational motion (z-translation)
            graph heave_;


        // Prescribed rotation members

            // - Graph for roll rotational velocity (x-rotation)
            graph roll_;

            // - Graph for pitch rotational velocity (y-rotation)
            graph pitch_;

            // - Graph for yaw rotational velocity (z-rotation)
            graph yaw_;

            //- Is the rotational velocity given in degrees/sec
            Switch inDegrees_;


        //- Member functions

            //- Return translational velocity vector
            vector translationalVelocity() const;

            // -Return rotational velocity vector
            vector rotationalVelocity() const;

            //- Calculate transformation septernion
            septernion calcTransformation() const;


public:

    //- Runtime type information
    TypeName("graphVelocity");


    // Constructors

        //- Construct from components
        graphVelocity
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );


    //- Destructor
    virtual ~graphVelocity()
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
