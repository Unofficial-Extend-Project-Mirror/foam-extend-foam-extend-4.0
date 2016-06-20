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

#include "linearOscillation.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(linearOscillation, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        linearOscillation,
        dictionary
    );
};
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector
Foam::solidBodyMotionFunctions::linearOscillation::calcPosition
(
    const scalar t
) const
{
    return amplitude_*sin(2*pi*t/period_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::linearOscillation::linearOscillation
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    amplitude_(SBMFCoeffs_.lookup("amplitude")),
    period_(readScalar(SBMFCoeffs_.lookup("period")))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::linearOscillation::~linearOscillation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::linearOscillation::transformation() const
{
    scalar t = time_.value();

    septernion TR(calcPosition(t), quaternion::I);

    Info<< "solidBodyMotionFunctions::linearOscillation::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


Foam::septernion
Foam::solidBodyMotionFunctions::linearOscillation::velocity() const
{
    scalar t = time_.value();
    scalar dt = time_.deltaT().value();

    septernion TV
    (
        (calcPosition(t) - calcPosition(t - dt))/dt,
        quaternion::zero
    );

    return TV;
}


bool Foam::solidBodyMotionFunctions::linearOscillation::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("amplitude") >> amplitude_;
    SBMFCoeffs_.lookup("amplitude") >> period_;

    return true;
}


// ************************************************************************* //
