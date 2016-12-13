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

#include "dampedCoulomb.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pairPotentials
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dampedCoulomb, 0);

addToRunTimeSelectionTable
(
    pairPotential,
    dampedCoulomb,
    dictionary
);

scalar dampedCoulomb::oneOverFourPiEps0 =
    1.0/(4.0 * mathematicalConstant::pi * 8.854187817e-12);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dampedCoulomb::dampedCoulomb
(
    const word& name,
    const dictionary& pairPotentialProperties
)
:
    pairPotential(name, pairPotentialProperties),
    dampedCoulombCoeffs_
    (
        pairPotentialProperties.subDict(typeName + "Coeffs")
    ),
    alpha_(readScalar(dampedCoulombCoeffs_.lookup("alpha")))
{
    setLookupTables();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar dampedCoulomb::unscaledEnergy(const scalar r) const
{
    return oneOverFourPiEps0*erfc(alpha_*r)/r;
}


bool dampedCoulomb::read(const dictionary& pairPotentialProperties)
{
    pairPotential::read(pairPotentialProperties);

    dampedCoulombCoeffs_ =
        pairPotentialProperties.subDict(typeName + "Coeffs");

    dampedCoulombCoeffs_.lookup("alpha") >> alpha_;

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pairPotentials
} // End namespace Foam

// ************************************************************************* //
