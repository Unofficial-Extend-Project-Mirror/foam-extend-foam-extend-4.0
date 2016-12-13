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
    Foam::solidBodyMotionFunctions::harmonicOscillation

Description
    Prescribed translational and rotational harmonic motion.
    The motion is of the following form: x = x_a*sin(omega*t + phi_x).
    Hence, velocities go with cosine instead of sine.

SourceFiles
    harmonicOscillation.C

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef harmonicOscillation_H
#define harmonicOscillation_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                         Class harmonicOscillation Declaration
\*---------------------------------------------------------------------------*/

class harmonicOscillation
:
    public solidBodyMotionFunction
{
    // Private data

        //- Centre of gravity
        point origin_;


        // Prescribed translation members

            //- Translational amplitudes
            vector transAmplitudes_;

            //- Translational angular frequencies
            vector transAngularFreq_;

            //- Translational phase angles
            vector transPhaseAngles_;


        // Prescribed rotation members

            //- Rotational amplitudes
            vector rotAmplitudes_;

            //- Rotational angular frequencies
            vector rotAngularFreq_;

            //- Rotational phase angles
            vector rotPhaseAngles_;

            //- Is the rotational velocity given in degrees/sec
            Switch inDegrees_;


    // Private Member Functions

        //- Disallow copy construct
        harmonicOscillation(const harmonicOscillation&);

        //- Disallow default bitwise assignment
        void operator=(const harmonicOscillation&);

        //- Calculate tranformation
        septernion calcTransformation(const scalar t) const;


public:

    //- Runtime type information
    TypeName("harmonicOscillation");


    // Constructors

        //- Construct from components
        harmonicOscillation
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );


    //- Destructor
    virtual ~harmonicOscillation()
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
