/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Foam::scalar Foam::solidBodyMotionFunctions::translation::rampFactor() const
{
    if (rampTime_ > 0)
    {
        // Calculate ramp-up factor
        const scalar t = time_.value();

        return sin(2*pi/(4*rampTime_)*Foam::min(rampTime_, t));
    }
    else
    {
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
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::translation::~translation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::translation::transformation() const
{
    scalar time = time_.value();

    septernion TR(rampFactor()*velocity_*time, quaternion::I);

    Info<< "solidBodyMotionFunctions::translation::transformation(): "
        << "Time = " << time << " velocity = " << rampFactor()*velocity_
        << " transformation = " << TR
        << endl;

    return TR;
}


Foam::septernion
Foam::solidBodyMotionFunctions::translation::velocity() const
{
    scalar time = time_.value();

    septernion TV(rampFactor()*velocity_, quaternion::zero);

    Info<< "solidBodyMotionFunctions::translation::transformation(): "
        << "Time = " << time << " velocity: " << TV << endl;

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
