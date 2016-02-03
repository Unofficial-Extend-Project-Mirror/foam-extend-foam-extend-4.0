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

#include "harmonicOscillation.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

defineTypeNameAndDebug(harmonicOscillation, 0);
addToRunTimeSelectionTable
(
    solidBodyMotionFunction,
    harmonicOscillation,
    dictionary
);

}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::harmonicOscillation::calcTransformation
(
    const scalar t
) const
{
    const vector translation =
    cmptMultiply
    (
        transAmplitudes_,
        vector
        (
            sin(transAngularFreq_.x()*t + transPhaseAngles_.x()),
            sin(transAngularFreq_.y()*t + transPhaseAngles_.y()),
            sin(transAngularFreq_.z()*t + transPhaseAngles_.z())
        )
    );


    const vector rotation =
    cmptMultiply
    (
        rotAmplitudes_,
        vector
        (
            sin(rotAngularFreq_.x()*t + rotPhaseAngles_.x()),
            sin(rotAngularFreq_.y()*t + rotPhaseAngles_.y()),
            sin(rotAngularFreq_.z()*t + rotPhaseAngles_.z())
        )
    );

    const   quaternion R(rotation.x(), rotation.y(), rotation.z());
    const   septernion TR
        (
            septernion(origin_ + translation)*R*septernion(-origin_)
        );

    return TR;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::harmonicOscillation::harmonicOscillation
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
Foam::solidBodyMotionFunctions::harmonicOscillation::transformation() const
{
   const scalar t = time_.value();

   const septernion TR = calcTransformation(t);

   Info<< "solidBodyMotionFunctions::harmonicOscillation::transformation(): "
   << "Time = " << t << " transformation: " << TR << endl;

   return TR;
}


Foam::septernion
Foam::solidBodyMotionFunctions::harmonicOscillation::velocity() const
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


bool Foam::solidBodyMotionFunctions::harmonicOscillation::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("origin") >> origin_;
    SBMFCoeffs_.lookup("translationalAmplitudes") >> transAmplitudes_;
    SBMFCoeffs_.lookup("translationalAngularFrequencies") >> transAngularFreq_;
    SBMFCoeffs_.lookup("translationalPhaseAngles") >> transPhaseAngles_;
    SBMFCoeffs_.lookup("rotationalAmplitudes") >> rotAmplitudes_;
    SBMFCoeffs_.lookup("rotationalAngularFrequencies") >> rotAngularFreq_;
    SBMFCoeffs_.lookup("rotationalPhaseAngles") >> rotPhaseAngles_;
    SBMFCoeffs_.lookup("inDegrees") >> inDegrees_;

    // Sanity check for negative frequencies
    if
    (
        transAngularFreq_.x() < -SMALL
     || transAngularFreq_.y() < -SMALL
     || transAngularFreq_.z() < -SMALL
     || rotAngularFreq_.x() < -SMALL
     || rotAngularFreq_.y() < -SMALL
     || rotAngularFreq_.z() < -SMALL
    )
    {
        FatalErrorIn
        (
           "harmonicOscillation::read"
             "(\n"
             "    const fvMesh& mesh\n"
             "    const dictionary& dict\n"
             ")"
        )   << "Negative angular frequency detected."
            << abort(FatalError);
    }

    // Convert to radians if necessary, printing data
    if (inDegrees_)
    {
        const scalar piBy180 = mathematicalConstant::pi/180.0;

        transPhaseAngles_ *= piBy180;
        rotAmplitudes_ *= piBy180;
        rotPhaseAngles_ *= piBy180;

        Info<< "Data in degrees converted:" << nl
            << "translationalPhaseAngles = " << transPhaseAngles_ << nl
            << "rotationalAmplitudes = " << rotAmplitudes_ << nl
            << "rotationalPhaseAngles = " << rotPhaseAngles_ << endl;
    }

    // Count prescribed rotations: only two harmonic rotations can be prescribed
    label nRot = 0;

    // Rotation
    mag(rotAmplitudes_.x()) > SMALL ? ++nRot : /* null expression */ 0 ;
    mag(rotAmplitudes_.y()) > SMALL ? ++nRot : /* null expression */ 0 ;
    mag(rotAmplitudes_.z()) > SMALL ? ++nRot : /* null expression */ 0 ;

    // Check if more than two rotations prescribed
    if (nRot > 2)
    {
        FatalErrorIn
        (
            "harmonicOscillation::read"
             "(\n"
             "    const fvMesh& mesh\n"
             "    const dictionary& dict\n"
             ")"
        )   << "Only two harmonic rotations can be prescribed. Prescribing "
            << "more than two yields unphysical results. "
            << abort(FatalError);
    }

    return true;
}


// ************************************************************************* //
