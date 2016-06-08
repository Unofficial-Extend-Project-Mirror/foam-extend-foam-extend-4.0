/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "timeVaryingCoulomb.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "frictionContactModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(timeVaryingCoulomb, 0);
    addToRunTimeSelectionTable(frictionLaw, timeVaryingCoulomb, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
timeVaryingCoulomb::timeVaryingCoulomb
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
timeVaryingCoulomb::timeVaryingCoulomb
(
    const timeVaryingCoulomb& fricLaw
)
:
    frictionLaw(fricLaw),
    frictionLawDict_(fricLaw.frictionLawDict_),
    frictionCoeffSeries_(fricLaw.frictionCoeffSeries_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

timeVaryingCoulomb::~timeVaryingCoulomb()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar timeVaryingCoulomb::slipTraction(const scalar pressure)
{
    return
        frictionCoeffSeries_
        (
            frictionModel().patch().boundaryMesh()
                .mesh().time().timeOutputValue()
        )*pressure;
}


void timeVaryingCoulomb::writeDict(Ostream& os) const
{
    word keyword("frictionLawDict");
    os.writeKeyword(keyword)
        << frictionCoeffSeries_;
}

// ************************************************************************* //

} // end of namespace foam
