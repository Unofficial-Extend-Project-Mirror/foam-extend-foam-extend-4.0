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

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "slicedVolFields.H"
#include "slicedSurfaceFields.H"
#include "immersedBoundaryFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedBoundaryFvPatch, 1);

    addToRunTimeSelectionTable(fvPatch, immersedBoundaryFvPatch, polyPatch);
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::debug::tolerancesSwitch
Foam::immersedBoundaryFvPatch::nonOrthogonalFactor_
(
    "immersedBoundaryNonOrthogonalFactor",
    0.1
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::immersedBoundaryFvPatch::makeCf(slicedSurfaceVectorField& Cf) const
{
    // Correct face centres by re-cutting and inserting the immersed patch
    Cf.reset(ibPolyPatch_.correctedFaceCentres());

    // Insert the patch data for the immersed boundary
    // Note: use the face centres from the stand-alone patch within the IB
    // HJ, 30/Nov/2017
    Cf.boundaryField()[index()].UList::operator=
    (
        vectorField::subField(ibPolyPatch_.ibPatch().faceCentres(), size())
    );
}


void Foam::immersedBoundaryFvPatch::makeSf(slicedSurfaceVectorField& Sf) const
{
    // Correct face centres by re-cutting and inserting the immersed patch
    Sf.reset(ibPolyPatch_.correctedFaceAreas());

    // Insert the patch data for the immersed boundary
    // Note: use the corrected face areas from immersed boundary instead of
    // the stand-alone patch areas within the IB
    // HJ, 30/Nov/2017
    Sf.boundaryField()[index()].UList::operator=
    (
        vectorField::subField(ibPolyPatch_.correctedIbPatchFaceAreas(), size())
    );
}


void Foam::immersedBoundaryFvPatch::makeC(slicedVolVectorField& C) const
{
    // Correct face centres by re-cutting and inserting the immersed patch
    C.reset
    (
        ibPolyPatch_.correctedCellCentres(),
        ibPolyPatch_.correctedFaceCentres()
    );

    // Insert the patch data for the immersed boundary
    // Note: use the face centres from the stand-alone patch within the IB
    // HJ, 30/Nov/2017
    C.boundaryField()[index()].UList::operator=
    (
        vectorField::subField(ibPolyPatch_.ibPatch().faceCentres(), size())
    );
}


void Foam::immersedBoundaryFvPatch::makeV(scalarField& V) const
{
    // Correct face centres by re-cutting and inserting the immersed patch
    V = ibPolyPatch_.correctedCellVolumes();
}


void Foam::immersedBoundaryFvPatch::updatePhi
(
    DimensionedField<scalar, volMesh>& V,
    DimensionedField<scalar, volMesh>& V0,
    surfaceScalarField& phi
) const
{
    // Correct face fluxes for cut area and insert the immersed patch fluxes

    const fvMesh& mesh = boundaryMesh().mesh();

    scalar deltaT = mesh.time().deltaT().value();
    scalar rDeltaT = 1.0/deltaT;

    // Multiply the raw mesh motion flux with the masking function
    scalarField ibAreaRatio =
        mag(ibPolyPatch_.correctedFaceAreas())/mag(mesh.faceAreas());

    // Scaling of internal mesh flux field should be done only for the current
    // ib patch to avoid scaling multiple times in case of multiple Ib patches
    // present. (IG 3/Dec/2018)

    // Scale internalField
    scalarField& phiIn = phi.internalField();

    const labelList& deadFaces = ibPolyPatch_.deadFaces();
    forAll (deadFaces, dfI)
    {
        const label faceI = deadFaces[dfI];
        if (mesh.isInternalFace(faceI))
        {
            phiIn[faceI] = scalar(0);
        }
    }
 
    const labelList& cutFaces = ibPolyPatch_.ibFaces();
    forAll (cutFaces, cfI)
    {
        const label faceI = cutFaces[cfI];
        if (mesh.isInternalFace(faceI))
        {
            phiIn[faceI] *= ibAreaRatio[faceI];
        }
    }

    // Scale all other patches
    forAll (mesh.boundary(), patchI)
    {
        if (!isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
        {
            phi.boundaryField()[patchI] *=
                mesh.boundary()[patchI].patchSlice(ibAreaRatio);
        }
    }

    // Immersed boundary patch
    // Calculate the mesh motion flux from the old and new coordinate of
    // triangular face centres and the time step dotted with the new face area
    phi.boundaryField()[index()] =
    (
        ibPolyPatch_.motionDistance()
      & ibPolyPatch_.correctedIbPatchFaceAreas()
    )*rDeltaT;

    // Check and adjust the immersed boundary space conservation law
    // The mesh motion fluxes come from the actual mesh motion or the motion
    // of the immersed boundary
    // The new cell volumes come from the current mesh configuration
    // The space conservation law will be satisfied by adjusting either
    // the old or the new cell volume.  HJ, 15/Dec/2017

    // First sum up all the fluxes
    scalarField divPhi(mesh.nCells(), 0);

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    forAll (owner, faceI)
    {
        divPhi[owner[faceI]] += phiIn[faceI];
        divPhi[neighbour[faceI]] -= phiIn[faceI];
    }

    // Add the mesh motion fluxes from all patches including immersed boundary
    forAll (mesh.boundary(), patchI)
    {
        const unallocLabelList& pFaceCells =
            mesh.boundary()[patchI].faceCells();

        const scalarField& pssf = phi.boundaryField()[patchI];

        // Check for size since uninitialised ib patches can have zero size at
        // this point (IG 7/Nov/2018)
        if (pssf.size() > 0)
        {
            forAll (pFaceCells, faceI)
            {
                divPhi[pFaceCells[faceI]] += pssf[faceI];
            }
        }
    }

    // Use corrected cell volume
    scalarField& newVols = V.field();
    scalarField& oldVols = V0.field();

    // Multiply by the time-step size and add new volume
    scalarField magDivPhi = mag((newVols - oldVols)*rDeltaT - divPhi);

    // Note:
    // The immersed boundary is now in the new position.  Therefore, some
    // cells that were cut are no longer in the contact with the IB, meaning
    // that ALL cells need to be checked and corrected
    // HJ, 22/Dec/2017
    forAll (magDivPhi, cellI)
    {
        // if (magDivPhi[cellI] > SMALL)
        if (magDivPhi[cellI] > 1e-40)
        {
            // Attempt to correct via old volume
            scalar corrOldVol = newVols[cellI] - divPhi[cellI]*deltaT;

            // Pout<< "Flux maneouvre for cell " << cellI << ": "
            //     << " error: " << magDivPhi[cellI]
            //     << " V: " << newVols[cellI]
            //     << " V0: " << oldVols[cellI]
            //     << " divPhi: " << divPhi[cellI];

            if (corrOldVol < SMALL)
            {
                // Update new volume because old volume cannot carry
                // the correction
                newVols[cellI] = oldVols[cellI] + divPhi[cellI]*deltaT;
            }
            else
            {
                oldVols[cellI] = corrOldVol;
            }

            // scalar corrDivMeshPhi =
            //     mag((newVols[cellI] - oldVols[cellI]) - divPhi[cellI]*deltaT);
            // Pout<< " Corrected: " << corrDivMeshPhi << endl;
        }
    }
}


void Foam::immersedBoundaryFvPatch::makeDeltaCoeffs
(
    fvsPatchScalarField& dc
) const
{
    const vectorField d = delta();

    dc = 1.0/max((nf() & d), 0.05*mag(d));
}


void Foam::immersedBoundaryFvPatch::makeCorrVecs(fvsPatchVectorField& cv) const
{
    // Set patch non-orthogonality correction to zero
//    cv = vector::zero;

    vectorField& cvIn = const_cast<vectorField&>(cv.internalField());

    const fvMesh& mesh = boundaryMesh().mesh();

    // Get face addressing
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    // Calculate volume fraction for all cells
    const scalarField gamma = mesh.V().field()/mesh.cellVolumes();

    // Visit all internal faces.  If a corrected volume fraction is smaller
    // than a threshold, reset non-orthogonality for the face
    forAll (neighbour, faceI)
    {
        if
        (
            gamma[owner[faceI]] < nonOrthogonalFactor_()
         || gamma[neighbour[faceI]] < nonOrthogonalFactor_()
        )
        {
            // Thin live cut.  Reset correction vectors
            cvIn[faceI] = vector::zero;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryFvPatch::immersedBoundaryFvPatch
(
    const polyPatch& patch,
    const fvBoundaryMesh& bm
)
:
    fvPatch(patch, bm),
    ibPolyPatch_(refCast<const immersedBoundaryPolyPatch>(patch)),
    mesh_(bm.mesh())
{
    ibPolyPatch_.clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::immersedBoundaryFvPatch::size() const
{
    // Immersed boundary patch size equals to the number of intersected cells
    // HJ, 28/Nov/2017

    // Note: asking for patch size triggers the cutting which involves
    // parallel communication.  This should be avoided under read/write, ie
    // when the ibPolyPatch_ is not initialised.
    // Initialisation happens when the fvMesh is initialised, which should be
    // sufficient
    //  HJ, 12/Dec/2018
    // if (!ibPolyPatch_.active())
    // {
    //     return 0;
    // }

    return ibPolyPatch_.ibCells().size();
}


const Foam::unallocLabelList&
Foam::immersedBoundaryFvPatch::faceCells() const
{
    return ibPolyPatch_.ibCells();
}


Foam::tmp<Foam::vectorField> Foam::immersedBoundaryFvPatch::nf() const
{
    // The algorithm has been changed because basic IB patch information
    // (nf and delta) is used in assembly of derived information
    // (eg. deltaCoeffs) and circular dependency needs to be avoided.
    // nf and delta vectors shall be calculated directly from the intersected
    // patch.  HJ, 21/Mar/2019

    return ibPolyPatch_.ibPatch().faceNormals();
}


Foam::tmp<Foam::vectorField> Foam::immersedBoundaryFvPatch::delta() const
{
    // Not strictly needed: this is for debug only.  HJ, 5/Apr/2019
    return ibPolyPatch_.ibPatch().faceCentres() - Cn();
}


// ************************************************************************* //
