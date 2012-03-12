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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "overlapGgiFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(overlapGgiFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, overlapGgiFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Make patch weighting factors
void Foam::overlapGgiFvPatch::makeWeights(scalarField& w) const
{
    // Calculation of weighting factors is performed from the master
    // position, using reconstructed shadow cell centres
    // HJ, 2/Aug/2007
    if (overlapGgiPolyPatch_.master())
    {
        vectorField n = nf();

        // Note: mag in the dot-product.
        // For all valid meshes, the non-orthogonality will be less that
        // 90 deg and the dot-product will be positive.  For invalid
        // meshes (d & s <= 0), this will stabilise the calculation
        // but the result will be poor.  HJ, 24/Aug/2011
        scalarField nfc =
            mag(n & (overlapGgiPolyPatch_.reconFaceCellCentres() - Cf()));

        w = nfc/(mag(n & (Cf() - Cn())) + nfc);
    }
    else
    {
        // Pick up weights from the master side
        scalarField masterWeights(shadow().size());
        shadow().makeWeights(masterWeights);

        w = interpolate(1 - masterWeights);
    }
}


// Make patch face - neighbour cell distances
void Foam::overlapGgiFvPatch::makeDeltaCoeffs(scalarField& dc) const
{
    if (overlapGgiPolyPatch_.master())
    {
        // Stabilised form for bad meshes.  HJ, 24/Aug/2011
        vectorField d = delta();

        dc = 1.0/max(nf() & d, 0.05*mag(d));
    }
    else
    {
        scalarField masterDeltas(shadow().size());
        shadow().makeDeltaCoeffs(masterDeltas);
        dc = interpolate(masterDeltas);
    }
}


void Foam::overlapGgiFvPatch::makeCorrVecs(vectorField& cv) const
{
    // Non-orthogonality correction on a ggi interface
    // MB, 7/April/2009

    // Calculate correction vectors on coupled patches
    const scalarField& patchDeltaCoeffs = deltaCoeffs();

    vectorField patchDeltas = delta();
    vectorField n = nf();
    cv = n - patchDeltas*patchDeltaCoeffs;
}


const Foam::overlapGgiFvPatch& Foam::overlapGgiFvPatch::shadow() const
{
    const fvPatch& p = this->boundaryMesh()[overlapGgiPolyPatch_.shadowIndex()];
    return refCast<const overlapGgiFvPatch>(p);
}


// Return delta (P to N) vectors across coupled patch
Foam::tmp<Foam::vectorField> Foam::overlapGgiFvPatch::delta() const
{
    if (overlapGgiPolyPatch_.master())
    {
        return overlapGgiPolyPatch_.reconFaceCellCentres() - Cn();
    }
    else
    {
//         vectorField masterDelta = shadow().Cn()
//             - overlapGgiPolyPatch_.shadow().reconFaceCellCentres();

//         return interpolate(masterDelta);

        return interpolate
        (
            shadow().Cn() - overlapGgiPolyPatch_.shadow().reconFaceCellCentres()
        );
    }
}


Foam::tmp<Foam::labelField> Foam::overlapGgiFvPatch::interfaceInternalField
(
    const unallocLabelList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::overlapGgiFvPatch::transfer
(
    const Pstream::commsTypes,
    const unallocLabelList& interfaceData
) const
{
    notImplemented
    (
        "overlapGgiFvPatchField<Type>::"
        "transfer(const unallocLabelList& interfaceData) const"
    );

    return labelField::null();
}


Foam::tmp<Foam::labelField> Foam::overlapGgiFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes,
    const unallocLabelList& iF
) const
{
    return shadow().patchInternalField(iF);
}


// ************************************************************************* //
