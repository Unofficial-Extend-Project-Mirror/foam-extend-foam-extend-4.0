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

#include "regionWiseOversetAdjustPhi.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "oversetMesh.H"
#include "fvc.H"
#include "processorPolyPatch.H"
#include "emptyOversetFvPatchFields.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::regionWiseOversetAdjustPhi
(
    surfaceScalarField& phi,
    const volVectorField& U
)
{
    const fvMesh& mesh = phi.mesh();

    // If the mesh is moving, adjustment needs to be calculated on
    // relative fluxes.  HJ, 13/Feb/2009
    if (mesh.moving())
    {
        fvc::makeRelative(phi, U);
    }

    // Get overset mesh
    const oversetMesh& om = oversetMesh::New(mesh);

    // Get addressing to fringe faces
    const labelList& fringeFaces = om.fringeFaces();
    const boolList& fringeFaceFlips = om.fringeFaceFlips();

    // Get internal owner-neighbour addressing
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    // Get region split to identify separate mesh components
    const labelList& regionID = om.regionID();

    // Incoming and outgoing region fluxes
    scalarField regionIn(om.regions().size(), 0);
    scalarField regionOut(om.regions().size(), 0);

    // Incoming and outgoing region fluxes through the fringe
    scalarField regionFringeIn(om.regions().size(), 0);
    scalarField regionFringeOut(om.regions().size(), 0);

    // Get flux internal field (for fringe faces later on)
    scalarField& phiIn = phi.internalField();

    // STAGE 1
    // Loop through boundary fields and collect the incoming/outgoing fluxes
    forAll (phi.boundaryField(), patchI)
    {
        // Get necessary references
        const fvPatchVectorField& Up = U.boundaryField()[patchI];
        const fvsPatchScalarField& phip = phi.boundaryField()[patchI];
        const unallocLabelList& fc = Up.patch().faceCells();

        // All coupled and emptyOverset patches should not be taken into account
        if
        (
            !Up.coupled()
         && !isA<emptyOversetFvPatchVectorField>(Up)
        )
        {
            // Note: no need to distinguish between fixed value fluxes or
            // adjustable fluxes as only the fringe fluxes are going to be
            // scaled
            forAll (phip, i)
            {
                // Get current region index
                const label curRegion = regionID[fc[i]];

                if (phip[i] < 0.0)
                {
                    regionIn[curRegion] -= phip[i];
                }
                else
                {
                    regionOut[curRegion] += phip[i];
                }
            }
        }
    }

    // STAGE 2
    // Loop through fringe faces and collect incoming/outgoing fluxes

    // Note: fringe faces contain both internal and processor boundary
    // faces. Make sure processor faces are not counted twice.
    forAll (fringeFaces, ffI)
    {
        // Get face information
        const label curFace = fringeFaces[ffI];
        const bool curFlip = fringeFaceFlips[ffI];

        if (mesh.isInternalFace(curFace))
        {
            // Internal face

            // Get region index
            const label curRegion = regionID[owner[curFace]];

            // Check whether owner and neighbour belong to the same region
            if (curRegion != regionID[neighbour[curFace]])
            {
                FatalErrorIn
                (
                    "void Foam::regionWiseOversetAdjustPhi\n"
                    "(\n"
                    "    surfaceScalarField& phi,\n"
                    "    volVectorField& U\n"
                    ")"
                )   << "Region index different for owner and neighbour "
                    << "of face " << curFace
                    << abort(FatalError);
            }

            // Get current face flux
            const scalar& curPhi = phiIn[curFace];

            if (curFlip)
            {
                if (curPhi > 0.0)
                {
                    // Flux going into the region (out of the fringe).
                    // Note that positive sign is kept.
                    regionFringeIn[curRegion] += curPhi;
                }
                else
                {
                    // Flux coming out of the region (into the fringe).
                    // Note reverted sign.
                    regionFringeOut[curRegion] -= curPhi;
                }
            }
            else
            {
                if (curPhi > 0.0)
                {
                    // Flux going out of the region (into the fringe).
                    // Note that positive sign is kept.
                    regionFringeOut[curRegion] += curPhi;
                }
                else
                {
                    // Flux going into the region (out of the fringe).
                    // Note reverted sign.
                    regionFringeIn[curRegion] -= curPhi;
                }
            }
        }
        else
        {
            // Processor boundary fringe face
            // Find patch and face
            const label patchI = mesh.boundaryMesh().whichPatch(curFace);
            const label faceI = mesh.boundaryMesh()[patchI].whichFace(curFace);

            // Sanity check
            if (patchI < 0)
            {
                FatalErrorIn
                (
                    "void Foam::regionWiseOversetAdjustPhi\n"
                    "(\n"
                    "    surfaceScalarField& phi,\n"
                    "    volVectorField& U\n"
                    ")"
                )   << "Cannot find patch for fringe face" << curFace
                    << abort(FatalError);
            }

            // Only account for processor face from the owner processor side
            if (isA<processorPolyPatch>(mesh.boundaryMesh()[patchI]))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>
                    (
                        mesh.boundaryMesh()[patchI]
                    );

                if (procPatch.owner())
                {
                    // Processor patch, master side

                    // Get region index
                    const label curRegion =
                        regionID[mesh.boundary()[patchI].faceCells()[faceI]];

                    // Get reference to the current face flux
                    const scalar& curPhi = phi.boundaryField()[patchI][faceI];

                    if (curFlip)
                    {
	                if (curPhi > 0.0)
	                {
	                    // Flux going into the region (out of the fringe).
	                    // Note that positive sign is kept.
	                    regionFringeIn[curRegion] += curPhi;
	                }
	                else
	                {
	                    // Flux coming out of the region (into the fringe).
	                    // Note reverted sign.
	                    regionFringeOut[curRegion] -= curPhi;
	                }
	            }
	            else
	            {
	                if (curPhi > 0.0)
	                {
	                    // Flux going out of the region (into the fringe).
	                    // Note that positive sign is kept.
	                    regionFringeOut[curRegion] += curPhi;
	                }
	                else
	                {
	                    // Flux going into the region (out of the fringe).
	                    // Note reverted sign.
	                    regionFringeIn[curRegion] -= curPhi;
	                }
                    }
                }
            }
            else
            {
                FatalErrorIn
                (
                    "void Foam::regionWiseOversetAdjustPhi\n"
                    "(\n"
                    "    surfaceScalarField& phi,\n"
                    "    volVectorField& U\n"
                    ")"
                )   << "Patch for fringe face" << curFace
                    << " is not of processor type"
                    << abort(FatalError);
            }
        }
    }

    // Do a global reduce of fluxes
    reduce(regionIn, sumOp<scalarField>());
    reduce(regionOut, sumOp<scalarField>());
    reduce(regionFringeIn, sumOp<scalarField>());
    reduce(regionFringeOut, sumOp<scalarField>());

    // This method will inevitably fail when regionFringeOut is zero, which can
    // happen if we do not have fully enclosed mesh inside a background mesh
    // (i.e. when we have a mesh which consists of left and right regions
    // connected through a planar fringe layer)
    if (min(regionFringeOut) < SMALL)
    {
         WarningIn
         (
             "void Foam::regionWiseOversetAdjustPhi\n"
             "(\n"
             "    surfaceScalarField& phi,\n"
             "    volVectorField& U\n"
             ")"
         )   << "Outflow through fringe in a certain region is zero. "
             << "Skipping flux adjustment for: " << phi.name()
             << endl;

         return;
    }

    // Calculate region flux correction. Note: only fluxes going out of the
    // region through the fringe will be scaled.
    scalarField regionFringeFluxScale =
        (regionIn - regionOut + regionFringeIn)/regionFringeOut;

    if (oversetMesh::debug)
    {
        const scalarField totalRegionIn = regionIn + regionFringeIn;
        const scalarField totalRegionOut = regionOut + regionFringeOut;

        Info<< "Region fringe balance for " << phi.name()
            << ": in = " << totalRegionIn
            << ", out = " << totalRegionOut
            << ", fringeIn = " << regionFringeIn
            << ", fringeOut = " << regionFringeOut
            << ", balance = " << totalRegionOut - totalRegionIn
            << ", fluxScale = " << regionFringeFluxScale
            << endl;
    }

    // Go through all fringe faces on each region and balance the fluxes
    forAll (fringeFaces, ffI)
    {
        const label curFace = fringeFaces[ffI];
        const bool curFlip = fringeFaceFlips[ffI];

        if (mesh.isInternalFace(curFace))
        {
            // Internal face

            // Get region index
            const label curRegion = regionID[owner[curFace]];

            // Get reference to the flux for scaling
            scalar& curPhi = phiIn[curFace];

            // Scale fringe flux going out of the region
            if (curFlip && (curPhi < 0))
            {
                curPhi *= regionFringeFluxScale[curRegion];
            }
            else if (!curFlip && (curPhi > 0))
            {
                curPhi *= regionFringeFluxScale[curRegion];
            }
        }
        else
        {
            // Processor boundary fringe face
            // Find patch and face
            const label patchI = mesh.boundaryMesh().whichPatch(curFace);
            const label faceI = mesh.boundaryMesh()[patchI].whichFace(curFace);

            // We need to scale both owner and neighbour fluxes because they
            // represent the same flux

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    mesh.boundaryMesh()[patchI]
                );

            // Get reference to the flux for scaling
            scalar& curPhi = phi.boundaryField()[patchI][faceI];

            if (procPatch.owner()) // Owner side
            {
                // Get region index
                const label curRegion =
                    regionID[mesh.boundary()[patchI].faceCells()[faceI]];

                // Scale fringe flux going out of the region
                if (curFlip && (curPhi < 0))
                {
                    curPhi *= regionFringeFluxScale[curRegion];
                }
                else if (!curFlip && (curPhi > 0))
                {
                    curPhi *= regionFringeFluxScale[curRegion];
                }
            }
            else // Neighbouring processor side
            {
                // Get region index
                const label curRegion =
                    regionID[mesh.boundary()[patchI].faceCells()[faceI]];

                // Scale fringe flux going out of the region
                // Note: change in curPhi sign and curFlip
                if (!curFlip && (curPhi > 0))
                {
                    curPhi *= regionFringeFluxScale[curRegion];
                }
                else if (curFlip && (curPhi < 0))
                {
                    curPhi *= regionFringeFluxScale[curRegion];
                }
            }
        }
    }

    // If the mesh is moving, adjustment needs to be calculated on
    // relative fluxes.  Now reverting to absolute fluxes.  HJ, 13/Feb/2009
    if (mesh.moving())
    {
        fvc::makeAbsolute(phi, U);
    }
}


// ************************************************************************* //
