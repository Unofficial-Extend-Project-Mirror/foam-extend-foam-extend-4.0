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
    Foam::solidBodyMotionFunction::graphMotion

Description
    Prescribed translational and rotational motion given graphs for surge,
    sway, heave, roll, pitch and yaw. Hence, the motion (parts of trajectory) is
    given in graphs, not the velocity.

SourceFiles
    graphMotion.C

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef graphMotion_H
#define graphMotion_H

#include "solidBodyMotionFunction.H"
#include "graph.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                         Class graphMotion Declaration
\*---------------------------------------------------------------------------*/

class graphMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- Centre of gravity
        point origin_;


        // Prescribed translation members

            // - Graph for surge translational motion (x-translation)
            graph surge_;

            // - Graph for sway translational motion (y-translation)
            graph sway_;

            // - Graph for heave translational motion (z-translation)
            graph heave_;


        // Prescribed rotation members

            // - Graph for roll rotational motion (x-rotation)
            graph roll_;

            // - Graph for pitch rotational motion (y-rotation)
            graph pitch_;

            // - Graph for yaw rotational motion (z-rotation)
            graph yaw_;

            //- Is the rotational velocity given in degrees/sec
            Switch inDegrees_;


        //- Member functions

            //- Calculate transformation septernion
            septernion calcTransformation(const scalar t) const;


public:

    //- Runtime type information
    TypeName("graphMotion");


    // Constructors

        //- Construct from components
        graphMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );


    //- Destructor
    virtual ~graphMotion()
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
