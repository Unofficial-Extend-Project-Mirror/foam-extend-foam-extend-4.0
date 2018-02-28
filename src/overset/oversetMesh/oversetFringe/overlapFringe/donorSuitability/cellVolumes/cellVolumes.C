/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "cellVolumes.H"
#include "overlapFringe.H"
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
    const overlapFringe& overlapFringeAlgorithm,
    const dictionary& dict
)
:
    donorSuitability(overlapFringeAlgorithm, dict)
{
    // Get local donor suitability function using cell volumes
    const scalarField& localDsf = overlapFringeAlgorithm.mesh().V().field();

    // Combine donor suitability function data across processors for parallel
    // run
    this->combineDonorSuitabilityFunction(localDsf);
}


// ************************************************************************* //
