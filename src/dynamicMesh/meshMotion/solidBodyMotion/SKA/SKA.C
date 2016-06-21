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

#include "SKA.H"
#include "addToRunTimeSelectionTable.H"
#include "Tuple2.H"
#include "IFstream.H"
#include "interpolateXY.H"
#include "mathematicalConstants.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(SKA, 0);
    addToRunTimeSelectionTable(solidBodyMotionFunction, SKA, dictionary);
};
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::SKA::calcTransformation(const scalar t) const
{
    if (t < times_[0])
    {
        FatalErrorIn
        (
            "solidBodyMotionFunctions::SKA::transformation()"
        )   << "current time (" << t
            << ") is less than the minimum in the data table ("
            << times_[0] << ')'
            << exit(FatalError);
    }

    if (t > times_[times_.size()-1])
    {
        FatalErrorIn
        (
            "solidBodyMotionFunctions::SKA::transformation()"
        )   << "current time (" << t
            << ") is greater than the maximum in the data table ("
            << times_[times_.size()-1] << ')'
            << exit(FatalError);
    }

    translationRotationVectors TRV = interpolateXY
    (
        t,
        times_,
        values_
    );

    // Convert the rotational motion from deg to rad
    TRV[1] *= pi/180.0;

    quaternion R(TRV[1].x(), TRV[1].y(), TRV[1].z());
    septernion TR(septernion(CofG_ + TRV[0])*R*septernion(-CofG_));

    return TR;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::SKA::SKA
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

Foam::solidBodyMotionFunctions::SKA::~SKA()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion Foam::solidBodyMotionFunctions::SKA::transformation() const
{
    scalar t = time_.value();

    septernion TR = calcTransformation(t);

    Info<< "solidBodyMotionFunctions::SKA::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


Foam::septernion Foam::solidBodyMotionFunctions::SKA::velocity() const
{
    // Velocity is calculated as a difference
    scalar t = time_.value();
    scalar dt = time_.deltaT().value();

    const septernion velocity
    (
        (calcTransformation(t).t() - calcTransformation(t - dt).t())/dt,
        (calcTransformation(t).r()/calcTransformation(t - dt).r())/dt
    );

    return velocity;
}


bool Foam::solidBodyMotionFunctions::SKA::read(const dictionary& SBMFCoeffs)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    // If the timeDataFileName has changed read the file

    fileName newTimeDataFileName(SBMFCoeffs_.lookup("timeDataFileName"));

    if (newTimeDataFileName != timeDataFileName_)
    {
        timeDataFileName_ = newTimeDataFileName;

        IFstream dataStream(timeDataFileName_);

        if (dataStream.good())
        {
            List<Tuple2<scalar, translationRotationVectors> > timeValues
            (
                dataStream
            );

            times_.setSize(timeValues.size());
            values_.setSize(timeValues.size());

            forAll(timeValues, i)
            {
                times_[i] = timeValues[i].first();
                values_[i] = timeValues[i].second();
            }
        }
        else
        {
            FatalErrorIn
            (
                "solidBodyMotionFunctions::SKA::read(const dictionary&)"
            )   << "Cannot open time data file " << timeDataFileName_
                << exit(FatalError);
        }
    }

    SBMFCoeffs_.lookup("CofG") >> CofG_;

    return true;
}


// ************************************************************************* //
