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

#include "constantVelocity.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

defineTypeNameAndDebug(constantVelocity, 0);
addToRunTimeSelectionTable
(
    solidBodyMotionFunction,
    constantVelocity,
    dictionary
);

}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::constantVelocity::calcTransformation
(
    const scalar t
) const
{
    const vector translation = transVelocity_*(t - startMotionTime_);
    const vector rotation = rotVelocity_*(t - startMotionTime_);

    const quaternion R(rotation.x(), rotation.y(), rotation.z());
    const septernion TR
    (
        septernion(origin_ + translation)*R*septernion(-origin_)
    );

    return TR;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::constantVelocity::constantVelocity
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime)
{
    read(SBMFCoeffs);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::constantVelocity::transformation() const
{
   const scalar t = time_.value();

   const septernion TR = calcTransformation(t);

   Info<< "solidBodyMotionFunctions::constantVelocity::transformation(): "
   << "Time = " << t << " transformation: " << TR << endl;

   return TR;
}


Foam::septernion
Foam::solidBodyMotionFunctions::constantVelocity::velocity() const
{
    const scalar t = time_.value();
    const scalar dt = time_.deltaT().value();

    const septernion velocity
    (
        (calcTransformation(t).t() - calcTransformation(t - dt).t())/dt,
        (calcTransformation(t).r()/calcTransformation(t - dt).r())/dt
    );

    return velocity;
}


bool Foam::solidBodyMotionFunctions::constantVelocity::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("origin") >> origin_;
    SBMFCoeffs_.lookup("translationalVelocity") >> transVelocity_;
    SBMFCoeffs_.lookup("rotationalVelocity") >> rotVelocity_;
    SBMFCoeffs_.lookup("startMotionTime") >> startMotionTime_;
    SBMFCoeffs_.lookup("inDegrees") >> inDegrees_;

    // Convert to radians if necessary
    rotVelocity_ *= inDegrees_ ? mathematicalConstant::pi/180.0 : 1;

    return true;
}

// ************************************************************************* //
