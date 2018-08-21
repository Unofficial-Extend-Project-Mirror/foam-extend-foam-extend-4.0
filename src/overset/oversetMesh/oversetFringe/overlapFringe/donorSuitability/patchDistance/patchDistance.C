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

#include "patchDistance.H"
#include "oversetFringe.H"
#include "oversetRegion.H"
#include "patchWave.H"
#include "polyPatchID.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace donorSuitability
{

defineTypeNameAndDebug(patchDistance, 0);
addToRunTimeSelectionTable
(
    donorSuitability,
    patchDistance,
    dictionary
);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::donorSuitability::patchDistance::patchDistance
(
    const oversetFringe& oversetFringeAlgorithm,
    const dictionary& dict
)
:
    donorSuitability(oversetFringeAlgorithm, dict)
{
    // Get reference to fvMesh
    const fvMesh& mesh = oversetFringeAlgorithm.mesh();

    // Get distance patch names for master and donor regions
    wordList masterRegionPatchNames =
        coeffDict().lookup("masterRegionDistancePatches");
    wordList donorRegionPatchNames =
        coeffDict().lookup("donorRegionDistancePatches");

    // Insert patch IDs into hash sets
    labelHashSet masterPatchIDs(masterRegionPatchNames.size());
    labelHashSet donorPatchIDs(donorRegionPatchNames.size());

    // Master region distance patches
    forAll(masterRegionPatchNames, patchI)
    {
        polyPatchID pID(masterRegionPatchNames[patchI], mesh.boundaryMesh());

        if (pID.active())
        {
            masterPatchIDs.insert(pID.index());
        }
        else
        {
            FatalErrorIn
            (
                "donorSuitability::"
                "patchDistance::patchDistance()"
            )   << "Cannot find distance patch named: "
                << masterRegionPatchNames[patchI]
                << " for master region." << nl
                << "Available patch names: " << mesh.boundaryMesh().names()
                << abort(FatalError);
        }
    }

    // Donor regions distance patches
    forAll(donorRegionPatchNames, patchI)
    {
        polyPatchID pID(donorRegionPatchNames[patchI], mesh.boundaryMesh());

        if (pID.active())
        {
            donorPatchIDs.insert(pID.index());
        }
        else
        {
            FatalErrorIn
            (
                "donorSuitability::"
                "patchDistance::patchDistance()"
            )   << "Cannot find distance patch named: "
                << donorRegionPatchNames[patchI]
                << " for donor regions." << nl
                << "Available patch names: " << mesh.boundaryMesh().names()
                << abort(FatalError);
        }
    }

    // Calculate distance from specified patches and do not correct for accurate
    // near wall distance (=false paramater)
    patchWave masterDistance(mesh, masterPatchIDs, false);
    patchWave donorDistance(mesh, donorPatchIDs, false);

    // Combine both master and donor distances into a single field
    scalarField localDsf = masterDistance.distance();
    localDsf = min(localDsf, donorDistance.distance());

    // Combine donor suitability function data across processors for parallel
    // run
    this->combineDonorSuitabilityFunction(localDsf);
}


// ************************************************************************* //
