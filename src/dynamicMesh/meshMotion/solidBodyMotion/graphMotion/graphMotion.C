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

\*---------------------------------------------------------------------------*/

#include "graphMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "interpolateXY.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

defineTypeNameAndDebug(graphMotion, 0);
addToRunTimeSelectionTable(solidBodyMotionFunction, graphMotion, dictionary);

}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::graphMotion::calcTransformation
(
    const scalar t
) const
{
    const vector translation
    (
        interpolateXY(t, surge_.x(), surge_.y()),
        interpolateXY(t, sway_.x(), sway_.y()),
        interpolateXY(t, heave_.x(), heave_.y())
    );

    vector rotation
    (
        interpolateXY(t, roll_.x(), roll_.y()),
        interpolateXY(t, pitch_.x(), pitch_.y()),
        interpolateXY(t, yaw_.x(), yaw_.y())
    );

    if (inDegrees_)
    {
        const scalar piBy180 = mathematicalConstant::pi/180.0;

        rotation *= piBy180;
    }

    const   quaternion R(rotation.x(), rotation.y(), rotation.z());
    const   septernion TR
        (
            septernion(origin_ + translation)*R*septernion(-origin_)
        );

    return TR;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::graphMotion::graphMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
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
Foam::solidBodyMotionFunctions::graphMotion::transformation() const
{
   const scalar t = time_.value();
   const septernion TR = calcTransformation(t);

   Info<< "solidBodyMotionFunctions::graphMotion::transformation(): "
   << "Time = " << t << " transformation: " << TR << endl;

   return TR;
}


Foam::septernion
Foam::solidBodyMotionFunctions::graphMotion::velocity() const
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


bool Foam::solidBodyMotionFunctions::graphMotion::read
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
