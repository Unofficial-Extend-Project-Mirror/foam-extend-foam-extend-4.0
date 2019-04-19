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

#include "faceArea.H"
#include "oversetFringe.H"
#include "oversetRegion.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace donorSuitability
{

defineTypeNameAndDebug(faceArea, 0);
addToRunTimeSelectionTable
(
    donorSuitability,
    faceArea,
    dictionary
);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::donorSuitability::faceArea::faceArea
(
    const oversetFringe& oversetFringeAlgorithm,
    const dictionary& dict
)
:
    donorSuitability(oversetFringeAlgorithm, dict)
{
    // Get fvMesh reference
    const fvMesh& mesh = oversetFringeAlgorithm.mesh();

    // Get local donor suitability function using minium face area of a cell
    scalarField localDsf(mesh.nCells(), GREAT);

    // Get necessary mesh data
    const scalarField& magSfIn = mesh.magSf().internalField();
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    // Note: only internal faces of the mesh are considered, there's no need to
    // loop through boundary faces
    forAll(magSfIn, faceI)
    {
        const label& own = owner[faceI];
        const label& nei = neighbour[faceI];

        localDsf[own] = min(localDsf[own], magSfIn[faceI]);
        localDsf[nei] = min(localDsf[nei], magSfIn[faceI]);
    }

    // Combine donor suitability function data across processors for parallel
    // run
    this->combineDonorSuitabilityFunction(localDsf);
}


// ************************************************************************* //
