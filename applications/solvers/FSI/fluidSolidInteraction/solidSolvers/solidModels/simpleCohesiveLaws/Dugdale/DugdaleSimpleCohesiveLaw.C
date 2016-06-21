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

Description
    Dugdale simple cohesive law.

\*---------------------------------------------------------------------------*/

#include "DugdaleSimpleCohesiveLaw.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DugdaleSimpleCohesiveLaw, 0);
    addToRunTimeSelectionTable
    (
        simpleCohesiveLaw,
        DugdaleSimpleCohesiveLaw,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::DugdaleSimpleCohesiveLaw::DugdaleSimpleCohesiveLaw
(
    const word& cohesiveLawName,
    const dictionary& dict
)
:
    simpleCohesiveLaw(cohesiveLawName, dict),
    deltaC_(GIc()/sigmaMax())
{}


Foam::DugdaleSimpleCohesiveLaw::DugdaleSimpleCohesiveLaw
(
    const DugdaleSimpleCohesiveLaw& dcl
)
:
    simpleCohesiveLaw(dcl),
    deltaC_(dcl.deltaC_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::DugdaleSimpleCohesiveLaw::~DugdaleSimpleCohesiveLaw()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return current holding traction
Foam::scalar Foam::DugdaleSimpleCohesiveLaw::traction(scalar delta) const
{
    if (delta > deltaC().value())
    {
        return 0.0;
    }
    else if (delta < 0)
    {
        return sigmaMax().value();
    }

    return sigmaMax().value();
}

// ************************************************************************* //
