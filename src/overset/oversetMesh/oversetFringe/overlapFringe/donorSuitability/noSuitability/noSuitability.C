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

#include "noSuitability.H"
#include "oversetFringe.H"
#include "oversetRegion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace donorSuitability
{

defineTypeNameAndDebug(noSuitability, 0);
addToRunTimeSelectionTable(donorSuitability, noSuitability, dictionary);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::donorSuitability::noSuitability::noSuitability
(
    const oversetFringe& oversetFringeAlgorithm,
    const dictionary& dict
)
:
    donorSuitability(oversetFringeAlgorithm, dict)
{
    // Need to initialise donor suitability function
    const scalarField localDsf(oversetFringeAlgorithm.mesh().nCells(), 0);
    this->combineDonorSuitabilityFunction(localDsf);

    // Set threshold to SMALL such that all the pairs become suitable
    this->threshold() = SMALL;
}


// ************************************************************************* //
