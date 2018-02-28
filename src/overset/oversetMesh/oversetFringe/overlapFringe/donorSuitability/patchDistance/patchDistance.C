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

#include "patchDistance.H"
#include "overlapFringe.H"
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
    const overlapFringe& overlapFringeAlgorithm,
    const dictionary& dict
)
:
    donorSuitability(overlapFringeAlgorithm, dict)
{
    // Get reference to fvMesh
    const fvMesh& mesh = overlapFringeAlgorithm.mesh();

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
