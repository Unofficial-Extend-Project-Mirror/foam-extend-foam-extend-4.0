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

#include "cellVolumes.H"
#include "oversetFringe.H"
#include "oversetRegion.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace donorSuitability
{

defineTypeNameAndDebug(cellVolumes, 0);
addToRunTimeSelectionTable
(
    donorSuitability,
    cellVolumes,
    dictionary
);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::donorSuitability::cellVolumes::cellVolumes
(
    const oversetFringe& oversetFringeAlgorithm,
    const dictionary& dict
)
:
    donorSuitability(oversetFringeAlgorithm, dict)
{
    // Get local donor suitability function using cell volumes
    const scalarField& localDsf = oversetFringeAlgorithm.mesh().V().field();

    // Combine donor suitability function data across processors for parallel
    // run
    this->combineDonorSuitabilityFunction(localDsf);
}


// ************************************************************************* //
