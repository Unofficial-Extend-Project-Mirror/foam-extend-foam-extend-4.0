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

#include "coulombOrowanFriction.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coulombOrowanFriction, 0);
    addToRunTimeSelectionTable(frictionLaw, coulombOrowanFriction, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
coulombOrowanFriction::coulombOrowanFriction
(
    const word& name,
    const frictionContactModel& fricModel,
    const dictionary& dict
)
:
    frictionLaw(name, fricModel, dict),
    frictionLawDict_(dict.subDict("frictionLawDict")),
    frictionCoeff_(readScalar(frictionLawDict_.lookup("frictionCoeff"))),
    shearStressLimit_(readScalar(frictionLawDict_.lookup("shearStressLimit")))
{}


// Construct as a copy
coulombOrowanFriction::coulombOrowanFriction
(
    const coulombOrowanFriction& fricLaw
)
:
    frictionLaw(fricLaw),
    frictionLawDict_(fricLaw.frictionLawDict_),
    frictionCoeff_(fricLaw.frictionCoeff_),
    shearStressLimit_(fricLaw.shearStressLimit_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

coulombOrowanFriction::~coulombOrowanFriction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar coulombOrowanFriction::slipTraction(const scalar pressure)
{
  // limit the shear
  return min(frictionCoeff_*pressure, shearStressLimit_);
}

void coulombOrowanFriction::writeDict(Ostream& os) const
{
  word keyword("frictionLawDict");
  os.writeKeyword(keyword)
    << frictionLawDict_;
}

// ************************************************************************* //

} // end of namespace foam
