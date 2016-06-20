/*---------------------------------------------------------------------------* \
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

#include "sixDoFRigidBodyMotionRestraint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::sixDoFRigidBodyMotionRestraint, 0);

defineRunTimeSelectionTable(Foam::sixDoFRigidBodyMotionRestraint, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraint::sixDoFRigidBodyMotionRestraint
(
    const dictionary& sDoFRBMRDict
)
:
    sDoFRBMRCoeffs_
    (
        sDoFRBMRDict.subDict
        (
            word(sDoFRBMRDict.lookup("sixDoFRigidBodyMotionRestraint"))
          + "Coeffs"
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraint::~sixDoFRigidBodyMotionRestraint()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::sixDoFRigidBodyMotionRestraint::read
(
    const dictionary& sDoFRBMRDict
)
{
    sDoFRBMRCoeffs_ = sDoFRBMRDict.subDict(type() + "Coeffs");

    return true;
}


// ************************************************************************* //
