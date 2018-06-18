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

#include "angularDamper.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(angularDamper, 0);
    addToRunTimeSelectionTable
    (
        rotationalRestraint,
        angularDamper,
        word
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::angularDamper::angularDamper
(
    const word& name,
    const dictionary& dict,
    const sixDOFODE& sixDOF
)
:
    rotationalRestraint(name, dict, sixDOF),
    angDampingCoeffs_(dict.lookup("angularDamping")),
    inGlobal_(dict.lookup("inGlobalCoordinateSystem"))
{}


Foam::autoPtr<Foam::rotationalRestraint>
Foam::angularDamper::clone() const
{
    return autoPtr<rotationalRestraint>
    (
        new angularDamper(*this)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::angularDamper::~angularDamper()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::angularDamper::restrainingMoment
(
    const scalar,
    const tensor& toRelative,
    const vector& omega
) const
{
    vector rm;

    if (inGlobal_)
    {
        // Restraint given in global (inertial) coordinate system, transform it
        // to local
        rm = (toRelative & angDampingCoeffs_) & omega;
    }
    else
    {
        // Restraint already in local (body) coordinate system
        rm = angDampingCoeffs_ & omega;
    }

    return -rm;
}


void Foam::angularDamper::write(Ostream& os) const
{
    os.writeKeyword("type") << tab << type()
        << token::END_STATEMENT << nl << nl;

    os.writeKeyword("angularDamping") << tab << angDampingCoeffs_
       << token::END_STATEMENT << nl;
    os.writeKeyword("inGlobalCoordinateSystem") << tab << inGlobal_
       << token::END_STATEMENT << endl;
}


// ************************************************************************* //
