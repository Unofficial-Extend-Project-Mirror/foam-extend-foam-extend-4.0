/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "periodicOscillation.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(periodicOscillation, 0);
    addToRunTimeSelectionTable
    (
        translationalConstraint,
        periodicOscillation,
        word
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::periodicOscillation::periodicOscillation
(
    const word& name,
    const dictionary& dict,
    const sixDOFODE& sixDOF
)
:
    translationalConstraint(name, dict, sixDOF),
    dir_(dict.lookup("direction")),
    a_(readScalar(dict.lookup("motionAmplitude"))),
    period_(readScalar(dict.lookup("period"))),
    omega_(mathematicalConstant::twoPi/period_),
    phi_(readScalar(dict.lookup("phaseShift"))*mathematicalConstant::pi/180.0)
{
    // Rescale direction
    if (mag(dir_) < SMALL)
    {
        FatalErrorIn
        (
            "Foam::periodicOscillation::periodicOscillation"
            "\n("
            "\n    const word& name,"
            "\n    const dictionary& dict"
            "\n)"
        )   << "Zero direction specified. This is not allowed."
            << exit(FatalError);
    }
    else
    {
        dir_ /= mag(dir_);
    }
}


Foam::autoPtr<Foam::translationalConstraint>
Foam::periodicOscillation::clone() const
{
    return autoPtr<translationalConstraint>
    (
        new periodicOscillation(*this)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::periodicOscillation::~periodicOscillation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::periodicOscillation::matrixContribution
(
    const scalar,
    const tensor&,
    const vector&,
    const vector&
) const
{
    return dir_;
}


Foam::scalar Foam::periodicOscillation::sourceContribution
(
    const scalar t,
    const tensor&,
    const vector&,
    const vector&
) const
{
    return -a_*sqr(omega_)*(sin(omega_*t + phi_));
}


void Foam::periodicOscillation::stabilise
(
    const scalar t,
    vector& x,
    vector& u
) const
{
    // Set the displacement according to periodic oscillation

    // First subtract calculated displacement...
    x -= (x & dir_)*dir_;

    // ... then add the correct displacement
    x += dir_*a_*sin(omega_*t + phi_);


    // Set the velocity according to periodic oscillation

    // First subract calculated velocity...
    u -= (u & dir_)*dir_;

    // ... then add the correct velocity
    u += dir_*a_*omega_*cos(omega_*t + phi_);
}


void Foam::periodicOscillation::write(Ostream& os) const
{
    os.writeKeyword("type") << tab << type()
        << token::END_STATEMENT << nl << nl;

    os.writeKeyword("direction") << tab << dir_
       << token::END_STATEMENT << nl;
    os.writeKeyword("motionAmplitude") << tab << a_
       << token::END_STATEMENT << nl;
    os.writeKeyword("period") << tab << period_
      << token::END_STATEMENT << nl;
    os.writeKeyword("phaseShift") << tab
      << phi_*180.0/mathematicalConstant::pi << token::END_STATEMENT << endl;
}


// ************************************************************************* //
