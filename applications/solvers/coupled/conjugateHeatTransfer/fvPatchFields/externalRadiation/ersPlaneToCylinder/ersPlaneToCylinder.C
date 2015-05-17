/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "ersPlaneToCylinder.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(ersPlaneToCylinder, 0);
addToRunTimeSelectionTable(externalRadiationSource, ersPlaneToCylinder, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ersPlaneToCylinder::ersPlaneToCylinder
(
    const word& name,
    const dictionary& dict,
    const fvPatch& p
)
:
    ersViewFactor(name, dict),
    direction_(dict.lookup("direction"))
{
    scalarField cosBeta = (direction_ & p.Sf()/p.magSf());

    F() = 0.5 - 0.5*cosBeta;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ersPlaneToCylinder::write(Ostream& os) const
{
    ersViewFactor::write(os);
    os.writeKeyword("direction") << direction_ << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
