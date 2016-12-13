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

#include "restrainedHarmonicSpring.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace tetherPotentials
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(restrainedHarmonicSpring, 0);

addToRunTimeSelectionTable
(
    tetherPotential,
    restrainedHarmonicSpring,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

restrainedHarmonicSpring::restrainedHarmonicSpring
(
    const word& name,
    const dictionary& tetherPotentialProperties
)
:
    tetherPotential(name, tetherPotentialProperties),
    restrainedHarmonicSpringCoeffs_
    (
        tetherPotentialProperties.subDict(typeName + "Coeffs")
    ),
    springConstant_
    (
        readScalar(restrainedHarmonicSpringCoeffs_.lookup("springConstant"))
    ),
    rR_
    (
        readScalar(restrainedHarmonicSpringCoeffs_.lookup("rR"))
    )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar restrainedHarmonicSpring::energy(const vector r) const
{
    scalar magR = mag(r);

    if (magR < rR_)
    {
        return 0.5 * springConstant_ * magSqr(r);
    }
    else
    {
        return 0.5 * springConstant_ * rR_ * rR_
          + springConstant_ * rR_ * (magR - rR_);
    }
}


vector restrainedHarmonicSpring::force(const vector r) const
{
    scalar magR = mag(r);

    if (magR < rR_)
    {
        return -springConstant_ * r;
    }
    else
    {
        return -springConstant_ * rR_ * r / magR;
    }
}


bool restrainedHarmonicSpring::read
(
    const dictionary& tetherPotentialProperties
)
{
    tetherPotential::read(tetherPotentialProperties);

    restrainedHarmonicSpringCoeffs_ =
        tetherPotentialProperties.subDict(typeName + "Coeffs");

    restrainedHarmonicSpringCoeffs_.lookup("springConstant") >> springConstant_;
    restrainedHarmonicSpringCoeffs_.lookup("rR") >> rR_;

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace tetherPotentials
} // End namespace Foam

// ************************************************************************* //
