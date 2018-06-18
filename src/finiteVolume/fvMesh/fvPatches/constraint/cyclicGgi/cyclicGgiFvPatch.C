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

Author
    Martin Beaudoin, Hydro-Quebec, (2008)

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "cyclicGgiFvPatch.H"
#include "fvMesh.H"
#include "fvPatchFields.H"
#include "fvPatchFields.H"
#include "fvsPatchFields.H"
#include "slicedSurfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicGgiFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, cyclicGgiFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Make mesh cell centres.  Moved from fvMeshGeometry
void Foam::cyclicGgiFvPatch::makeC(slicedSurfaceVectorField& C) const
{
    C.boundaryField()[index()].UList<vector>::operator=
    (
        patchSlice(cyclicGgiPolyPatch_.boundaryMesh().mesh().faceCentres())
    );
}


// Make patch weighting factors
void Foam::cyclicGgiFvPatch::makeWeights(fvsPatchScalarField& w) const
{
    // Calculation of weighting factors is performed from the master
    // position, using reconstructed shadow cell centres
    // HJ, 2/Aug/2007
    if (cyclicGgiPolyPatch_.master())
    {
        vectorField n = nf();

        // Note: mag in the dot-product.
        // For all valid meshes, the non-orthogonality will be less than
        // 90 deg and the dot-product will be positive.  For invalid
        // meshes (d & s <= 0), this will stabilise the calculation
        // but the result will be poor.  HJ, 24/Aug/2011
        scalarField nfc =
            mag
            (
                n & (cyclicGgiPolyPatch_.reconFaceCellCentres() - Cf())
            );

        w = nfc/(mag(n & (Cf() - Cn())) + nfc + SMALL);
    }
    else
    {
        // Slave side. Interpolate the master side weights, scale them for
        // partially covered faces and set weights for fully uncovered faces if
        // the bridge overlap is switched on. VV, 15/Feb/2018.

        // Pick up weights from the master side
        fvsPatchScalarField masterWeights
        (
            shadow(),
            w.dimensionedInternalField()
        );

        shadow().makeWeights(masterWeights);

        // Interpolate master weights to this side
        w = interpolate(masterWeights);

        if (bridgeOverlap())
        {
            // Weights for fully uncovered faces
            const scalarField uncoveredWeights(w.size(), 0.5);

            // Set weights for uncovered faces
            setUncoveredFaces(uncoveredWeights, w);

            // Scale partially overlapping faces
            scalePartialFaces(w);
        }

        // Finally construct these weights as 1 - master weights
        w = 1 - w;
    }
}


// Make patch face - neighbour cell distances
void Foam::cyclicGgiFvPatch::makeDeltaCoeffs(fvsPatchScalarField& dc) const
{
    if (cyclicGgiPolyPatch_.master())
    {
        // Master side. No need to scale partially uncovered or set fully
        // uncovered faces since delta already takes it into account.
        // VV, 25/Feb/2018.

        // Stabilised form for bad meshes.  HJ, 24/Aug/2011
        const vectorField d = delta();

        dc = 1.0/max(nf() & d, 0.05*mag(d));
    }
    else
    {
        // Slave side. Interpolate the master side, scale it for partially
        // covered faces and set deltaCoeffs for fully uncovered faces if the
        // bridge overlap is switched on. VV, 15/Feb/2018.

        fvsPatchScalarField masterDeltas
        (
            shadow(),
            dc.dimensionedInternalField()
        );

        shadow().makeDeltaCoeffs(masterDeltas);

        dc = interpolate(masterDeltas);

        if (bridgeOverlap())
        {
            // Delta coeffs for fully uncovered faces obtained from deltas on
            // this side
            const vectorField d = delta();
            const scalarField uncoveredDeltaCoeffs =
                1.0/max(nf() & d, 0.05*mag(d));

            // Set delta coeffs for uncovered faces
            setUncoveredFaces(uncoveredDeltaCoeffs, dc);

            // Scale partially overlapping faces
            scalePartialFaces(dc);
        }
    }
}


const Foam::cyclicGgiFvPatch& Foam::cyclicGgiFvPatch::shadow() const
{
    const fvPatch& p = this->boundaryMesh()[cyclicGgiPolyPatch_.shadowIndex()];
    return refCast<const cyclicGgiFvPatch>(p);
}


// Return delta (P to N) vectors across coupled patch
Foam::tmp<Foam::vectorField> Foam::cyclicGgiFvPatch::delta() const
{
    if (cyclicGgiPolyPatch_.master())
    {
        // Master side. Note: scaling partially covered faces and setting deltas
        // to fully uncovered faces correctly taken into account in
        // reconFaceCellCentres function. VV, 15/Feb/2018.

        tmp<vectorField> tdelta =
            cyclicGgiPolyPatch_.reconFaceCellCentres() - Cn();

        return tdelta;
    }
    else
    {
        // Slave side. Interpolate the master side, scale it for partially
        // covered faces and set deltas for fully uncovered faces if the bridge
        // overlap is switched on. VV, 15/Feb/2018.

        tmp<vectorField> tdelta = interpolate
        (
            shadow().Cn() - cyclicGgiPolyPatch_.shadow().reconFaceCellCentres()
        );

        if (bridgeOverlap())
        {
            // Deltas for fully uncovered faces
            const vectorField uncoveredDeltas(2.0*fvPatch::delta());

            // Set deltas for fully uncovered faces
            setUncoveredFaces(uncoveredDeltas, tdelta());

            // Scale for partially covered faces
            scalePartialFaces(tdelta());
        }

        return tdelta;
    }
}


// ************************************************************************* //
