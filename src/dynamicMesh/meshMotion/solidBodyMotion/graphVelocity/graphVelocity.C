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

#include "graphVelocity.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "interpolateXY.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

defineTypeNameAndDebug(graphVelocity, 0);
addToRunTimeSelectionTable(solidBodyMotionFunction, graphVelocity, dictionary);

}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector
Foam::solidBodyMotionFunctions::graphVelocity::translationalVelocity() const
{
    const scalar t = time_.value() - time_.deltaT().value()/2;

    return vector
           (
               interpolateXY(t, surge_.x(), surge_.y()),
               interpolateXY(t, sway_.x(), sway_.y()),
               interpolateXY(t, heave_.x(), heave_.y())
           );
}


Foam::vector
Foam::solidBodyMotionFunctions::graphVelocity::rotationalVelocity() const
{
    const scalar t = time_.value() - time_.deltaT().value()/2;

    scalar rollVelocity = interpolateXY(t, roll_.x(), roll_.y());
    scalar pitchVelocity = interpolateXY(t, pitch_.x(), pitch_.y());
    scalar yawVelocity = interpolateXY(t, yaw_.x(), yaw_.y());

    if (inDegrees_)
    {
        const scalar piBy180 = mathematicalConstant::pi/180.0;

        rollVelocity *= piBy180;
        pitchVelocity *= piBy180;
        yawVelocity *= piBy180;
    }

    return vector(rollVelocity, pitchVelocity, yawVelocity);
}


Foam::septernion
Foam::solidBodyMotionFunctions::graphVelocity::calcTransformation() const
{
    // Integrate velocity to get position
    if(localTimeIndex_ != time_.timeIndex())
    {
        // Set old translation and rotation to the previous state
        translationOld_ = translation_;
        rotationOld_ = rotation_;

        const scalar dt = time_.deltaT().value();

        translation_ += translationalVelocity()*dt;
        rotation_ += rotationalVelocity()*dt;

        localTimeIndex_ = time_.timeIndex();
    }

    const quaternion R(rotation_.x(), rotation_.y(), rotation_.z());
    const septernion TR
        (
            septernion(origin_ + translation_)*R*septernion(-origin_)
        );

    return TR;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::graphVelocity::graphVelocity
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    translation_(vector::zero),
    rotation_(vector::zero),
    translationOld_(vector::zero),
    rotationOld_(vector::zero),
    localTimeIndex_(-1),
    surge_
    (
        "surge",
        "t",
        "eta1Dot",
        IFstream
        (
            time_.path()/time_.caseConstant()/
            word(SBMFCoeffs_.lookup("surge"))
        )()
    ),
    sway_
    (
        "sway",
        "t",
        "eta2Dot",
        IFstream
        (
            time_.path()/time_.caseConstant()/
            word(SBMFCoeffs_.lookup("sway"))
        )()
    ),
    heave_
    (
        "heave",
        "t",
        "eta3Dot",
        IFstream
        (
            time_.path()/time_.caseConstant()/
            word(SBMFCoeffs_.lookup("heave"))
        )()
    ),
    roll_
    (
        "roll",
        "t",
        "eta4Dot",
        IFstream
        (
            time_.path()/time_.caseConstant()/
            word(SBMFCoeffs_.lookup("roll"))
        )()
    ),
    pitch_
    (
        "pitch",
        "t",
        "eta5Dot",
        IFstream
        (
            time_.path()/time_.caseConstant()/
            word(SBMFCoeffs_.lookup("pitch"))
        )()
    ),
    yaw_
    (
        "yaw",
        "t",
        "eta6Dot",
        IFstream
        (
            time_.path()/time_.caseConstant()/
            word(SBMFCoeffs_.lookup("yaw"))
        )()
    )
{
    read(SBMFCoeffs);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::graphVelocity::transformation() const
{
   const septernion TR = calcTransformation();

   Info<< "solidBodyMotionFunctions::graphVelocity::transformation(): "
   << "Time = " << time_.value() << " transformation: " << TR << endl;

   return TR;
}


Foam::septernion
Foam::solidBodyMotionFunctions::graphVelocity::velocity() const
{
    const scalar dt = time_.deltaT().value();

    const septernion TR = calcTransformation();

    const quaternion ROld(rotationOld_.x(), rotationOld_.y(), rotationOld_.z());
    const septernion TROld
        (
            septernion(origin_ + translationOld_)*ROld*septernion(-origin_)
        );

    return septernion
    (
        (TR.t() - TROld.t())/dt,
        (TR.r()/TROld.r())/dt
    );
}


bool Foam::solidBodyMotionFunctions::graphVelocity::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("origin") >> origin_;
    SBMFCoeffs_.lookup("inDegrees") >> inDegrees_;

    word surge = SBMFCoeffs_.lookup("surge");
    word sway = SBMFCoeffs_.lookup("sway");
    word heave = SBMFCoeffs_.lookup("heave");
    word roll = SBMFCoeffs_.lookup("roll");
    word pitch = SBMFCoeffs_.lookup("pitch");
    word yaw = SBMFCoeffs_.lookup("yaw");

    surge_ = graph
    (
        "surge",
        "t",
        "eta1Dot",
        IFstream
        (
            time_.path()/time_.caseConstant()/
            surge
        )()
    );
    sway_ = graph
    (
        "sway",
        "t",
        "eta2Dot",
        IFstream
        (
            time_.path()/time_.caseConstant()/
            sway
        )()
    );
    heave_ = graph
    (
        "heave",
        "t",
        "eta3Dot",
        IFstream
        (
            time_.path()/time_.caseConstant()/
            heave
        )()
    );
    roll_ = graph
    (
        "roll",
        "t",
        "eta4Dot",
        IFstream
        (
            time_.path()/time_.caseConstant()/
            roll
        )()
    );
    pitch_ = graph
    (
        "pitch",
        "t",
        "eta5Dot",
        IFstream
        (
            time_.path()/time_.caseConstant()/
            pitch
        )()
    );
    yaw_ = graph
    (
        "yaw",
        "t",
        "eta6Dot",
        IFstream
        (
            time_.path()/time_.caseConstant()/
            yaw
        )()
    );

    return true;
}

// ************************************************************************* //
