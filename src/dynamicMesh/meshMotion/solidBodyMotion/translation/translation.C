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

\*---------------------------------------------------------------------------*/

#include "translation.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(translation, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        translation,
        dictionary
    );
};
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar
Foam::solidBodyMotionFunctions::translation::rampFactor() const
{
    const scalar t = time_.value();

    if (t < rampTime_)
    {
        // Ramping region
        return sin(pi/(2*rampTime_)*t);
    }
    else
    {
        // Past ramping region
        return 1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::translation::translation
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    velocity_(SBMFCoeffs_.lookup("velocity")),
    rampTime_(readScalar(SBMFCoeffs_.lookup("rampTime")))
{
    if (rampTime_ < 0)
    {
        FatalIOErrorIn
        (
            "solidBodyMotionFunctions::translation::translation",
            SBMFCoeffs_
        )   << "Negative rampTime not allowed."
            << abort(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::translation::~translation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::translation::transformation() const
{
    const scalar t = time_.value();

    septernion TR;

    if (t < rampTime_)
    {
        // Ramping region
        // Account for ramping using analytical integration from ramped velocity
        // distribution in this region
        TR = septernion
        (
            velocity_*2*rampTime_/pi*(1 - cos(pi/(2*rampTime_)*t)),
            quaternion::I
        );
    }
    else
    {
        // Past ramping region
        TR = septernion
        (
            velocity_*
            (
                // Displacement during the ramping region
                2*rampTime_/pi
                // Displacement during constant velocity after ramping region
              + (t - rampTime_)
            ),
            quaternion::I
        );
    }

    Info<< "solidBodyMotionFunctions::translation::transformation(): "
        << "Time = " << t << " velocity = " << rampFactor()*velocity_
        << " transformation = " << TR
        << endl;

    return TR;
}


Foam::septernion
Foam::solidBodyMotionFunctions::translation::velocity() const
{
    septernion TV
    (
        rampFactor()*velocity_,
        quaternion::I/time_.deltaT().value()
    );

    Info<< "solidBodyMotionFunctions::translation::transformation(): "
        << "Time = " << time_.value() << " velocity: " << TV << endl;

    return TV;
}


bool Foam::solidBodyMotionFunctions::translation::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    rampTime_ = readScalar(SBMFCoeffs_.lookup("rampTime"));

    return true;
}


// ************************************************************************* //
