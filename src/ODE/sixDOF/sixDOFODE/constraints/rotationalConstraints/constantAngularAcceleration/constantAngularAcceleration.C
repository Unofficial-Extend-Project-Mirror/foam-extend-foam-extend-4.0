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

#include "constantAngularAcceleration.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantAngularAcceleration, 0);
    addToRunTimeSelectionTable
    (
        rotationalConstraint,
        constantAngularAcceleration,
        word
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantAngularAcceleration::constantAngularAcceleration
(
    const word& name,
    const dictionary& dict,
    const sixDOFODE& sixDOF
)
:
    rotationalConstraint(name, dict, sixDOF),
    dir_(dict.lookup("constraintDirection")),
    alpha_(readScalar(dict.lookup("angularAcceleration"))),
    inGlobal_(dict.lookup("inGlobalCoordinateSystem"))
{
    // Rescale direction
    if (mag(dir_) < SMALL)
    {
        FatalErrorIn
        (
            "Foam::constantTranslationalAcceleration::"
            "constantTranslationalAcceleration"
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


Foam::autoPtr<Foam::rotationalConstraint>
Foam::constantAngularAcceleration::clone() const
{
    return autoPtr<rotationalConstraint>
    (
        new constantAngularAcceleration(*this)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantAngularAcceleration::~constantAngularAcceleration()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::constantAngularAcceleration::matrixContribution
(
    const scalar,
    const tensor& toRelative,
    const vector&
) const
{
    vector mc;

    if (inGlobal_)
    {
        // Constraint is given in global (inertial) coordinate system, transform
        // it to local
        mc = toRelative & dir_;
    }
    else
    {
        // Constraint already in local (body) coordinate system
        mc = dir_;
    }

    return mc;
}


Foam::scalar Foam::constantAngularAcceleration::sourceContribution
(
    const scalar,
    const tensor& toRelative,
    const vector&
) const
{
    return alpha_;
}


void Foam::constantAngularAcceleration::write(Ostream& os) const
{
    os.writeKeyword("type") << tab << type()
       << token::END_STATEMENT << nl << nl;

    os.writeKeyword("constraintDirection") << tab << dir_
       << token::END_STATEMENT << nl;
    os.writeKeyword("angularAcceleration") << tab << alpha_
       << token::END_STATEMENT << nl;
    os.writeKeyword("inGlobalCoordinateSystem") << tab << inGlobal_
       << token::END_STATEMENT << endl;
}


// ************************************************************************* //
