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

#include "pressureVaryingCoulomb.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "frictionContactModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pressureVaryingCoulomb, 0);
    addToRunTimeSelectionTable(frictionLaw, pressureVaryingCoulomb, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
pressureVaryingCoulomb::pressureVaryingCoulomb
(
    const word& name,
    const frictionContactModel& fricModel,
    const dictionary& dict
)
:
    frictionLaw(name, fricModel, dict),
    frictionLawDict_(dict.subDict("frictionLawDict")),
    frictionCoeffSeries_(frictionLawDict_)
{}


// Construct as a copy
pressureVaryingCoulomb::pressureVaryingCoulomb
(
    const pressureVaryingCoulomb& fricLaw
)
:
    frictionLaw(fricLaw),
    frictionLawDict_(fricLaw.frictionLawDict_),
    frictionCoeffSeries_(fricLaw.frictionCoeffSeries_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pressureVaryingCoulomb::~pressureVaryingCoulomb()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar pressureVaryingCoulomb::slipTraction(const scalar pressure)
{
    return
        frictionCoeffSeries_
        (
            (pressure/1.0e9)
        )*pressure;
}


void pressureVaryingCoulomb::writeDict(Ostream& os) const
{
    word keyword("frictionLawDict");
    os.writeKeyword(keyword)
        << frictionCoeffSeries_;
}

// ************************************************************************* //

} // end of namespace foam
