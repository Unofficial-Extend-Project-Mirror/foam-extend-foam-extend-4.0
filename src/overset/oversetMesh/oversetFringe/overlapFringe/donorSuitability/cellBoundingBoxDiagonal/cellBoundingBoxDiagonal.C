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

#include "cellBoundingBoxDiagonal.H"
#include "overlapFringe.H"
#include "oversetRegion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace donorSuitability
{

defineTypeNameAndDebug(cellBoundingBoxDiagonal, 0);
addToRunTimeSelectionTable
(
    donorSuitability,
    cellBoundingBoxDiagonal,
    dictionary
);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::donorSuitability::cellBoundingBoxDiagonal::cellBoundingBoxDiagonal
(
    const overlapFringe& overlapFringeAlgorithm,
    const dictionary& dict
)
:
    donorSuitability(overlapFringeAlgorithm, dict)
{
    // Get reference to fvMesh
    const fvMesh& mesh = overlapFringeAlgorithm.mesh();

    // Get necessary mesh data
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();

    // Create local donor suitability function
    scalarField localDsf(mesh.nCells(), 0);

    // Loop through cells and calculate the bounding box diagonal of each cell
    forAll (cells, cellI)
    {
        const boundBox bb(cells[cellI].points(faces, points), false);
        localDsf[cellI] = bb.mag();
    }

    // Combine donor suitability function data across processors for parallel
    // run
    this->combineDonorSuitabilityFunction(localDsf);
}


// ************************************************************************* //
