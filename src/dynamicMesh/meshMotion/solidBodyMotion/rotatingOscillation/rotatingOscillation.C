/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "rotatingOscillation.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(rotatingOscillation, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        rotatingOscillation,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector
Foam::solidBodyMotionFunctions::rotatingOscillation::calcPosition
(
    const scalar t
) const
{
    return amplitude_*sin(2*pi*t/period_);
}


Foam::septernion
Foam::solidBodyMotionFunctions::rotatingOscillation::calcTransformation
(
    const scalar t
) const
{
    vector eulerAngles = amplitude_*sin(2*pi*t/period_);

    // Convert the rotational motion from deg to rad
    eulerAngles *= pi/180.0;

    const   quaternion R(eulerAngles.x(), eulerAngles.y(), eulerAngles.z());
    const   septernion TR(septernion(origin_)*R*septernion(-origin_));

    return TR;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::rotatingOscillation::
rotatingOscillation
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime)
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::rotatingOscillation::
~rotatingOscillation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::rotatingOscillation::
transformation() const
{
    scalar t = time_.value();

    const septernion TR = calcTransformation(t);

    Info<< "solidBodyMotionFunctions::rotatingOscillation::"
        << "transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


Foam::septernion
Foam::solidBodyMotionFunctions::rotatingOscillation::velocity() const
{
    scalar t = time_.value();
    scalar dt = time_.deltaT().value();

    const septernion velocity
    (
        (calcTransformation(t).t() - calcTransformation(t - dt).t())/dt,
        (calcTransformation(t).r()/calcTransformation(t - dt).r())/dt
    );

    return velocity;
}


bool Foam::solidBodyMotionFunctions::rotatingOscillation::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("origin") >> origin_;
    SBMFCoeffs_.lookup("amplitude") >> amplitude_;
    SBMFCoeffs_.lookup("period") >> period_;

    return true;
}


// ************************************************************************* //
