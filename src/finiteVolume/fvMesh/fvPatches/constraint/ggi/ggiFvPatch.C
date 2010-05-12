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

Description
    Generalized grid interface (GGI) patch, providing coupling
    between arbitrary patches which belong to the same fvMesh

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Contributor
    Martin Beaudoin, Hydro-Quebec, (2008)

\*---------------------------------------------------------------------------*/

#include "ggiFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ggiFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, ggiFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ggiFvPatch::~ggiFvPatch()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Make patch weighting factors
void Foam::ggiFvPatch::makeWeights(scalarField& w) const
{
    // Calculation of weighting factors is performed from the master
    // position, using reconstructed shadow cell centres
    // HJ, 2/Aug/2007
    if (ggiPolyPatch_.master())
    {
        vectorField n = nf();
        scalarField nfc = n & (ggiPolyPatch_.reconFaceCellCentres() - Cf());

        w = nfc/((n & (Cf() - Cn())) + nfc);

        if (bridgeOverlap())
        {
            // Set overlap weights to 0.5 and use mirrored neighbour field
            // for interpolation.  HJ, 21/Jan/2009
            bridge(scalarField(size(), 0.5), w);
        }
    }
    else
    {
        // Pick up weights from the master side
        scalarField masterWeights(shadow().size());
        shadow().makeWeights(masterWeights);

        w = interpolate(1 - masterWeights);

        if (bridgeOverlap())
        {
            // Set overlap weights to 0.5 and use mirrored neighbour field
            // for interpolation.  HJ, 21/Jan/2009
            bridge(scalarField(size(), 0.5), w);
        }
    }
}


// Make patch face - neighbour cell distances
void Foam::ggiFvPatch::makeDeltaCoeffs(scalarField& dc) const
{
    if (ggiPolyPatch_.master())
    {
        dc = 1.0/max(nf() & delta(), 0.05*mag(delta()));

        if (bridgeOverlap())
        {
            scalarField bridgeDeltas = nf() & fvPatch::delta();

            bridge(bridgeDeltas, dc);
        }
    }
    else
    {
        scalarField masterDeltas(shadow().size());
        shadow().makeDeltaCoeffs(masterDeltas);
        dc = interpolate(masterDeltas);

        if (bridgeOverlap())
        {
            scalarField bridgeDeltas = nf() & fvPatch::delta();

            bridge(bridgeDeltas, dc);
        }
    }
}


void Foam::ggiFvPatch::makeCorrVecs(vectorField& cv) const
{
    // Non-orthogonality correction on a ggi interface
    // MB, 7/April/2009

    // Calculate correction vectors on coupled patches
    const scalarField& patchDeltaCoeffs = deltaCoeffs();

    vectorField patchDeltas = delta();
    vectorField n = nf();
    cv = n - patchDeltas*patchDeltaCoeffs;
}


const Foam::ggiFvPatch& Foam::ggiFvPatch::shadow() const
{
    const fvPatch& p = this->boundaryMesh()[ggiPolyPatch_.shadowIndex()];

    return refCast<const ggiFvPatch>(p);
}


// Return delta (P to N) vectors across coupled patch
Foam::tmp<Foam::vectorField> Foam::ggiFvPatch::delta() const
{
    if (ggiPolyPatch_.master())
    {
        tmp<vectorField> tDelta = ggiPolyPatch_.reconFaceCellCentres() - Cn();

        if (bridgeOverlap())
        {
            vectorField bridgeDeltas = Cf() - Cn();

            bridge(bridgeDeltas, tDelta());
        }

        return tDelta;
    }
    else
    {
        tmp<vectorField> tDelta = interpolate
        (
            shadow().Cn() - ggiPolyPatch_.shadow().reconFaceCellCentres()
        );

        if (bridgeOverlap())
        {
            vectorField bridgeDeltas = Cf() - Cn();

            bridge(bridgeDeltas, tDelta());
        }

        return tDelta;
    }
}


bool Foam::ggiFvPatch::master() const
{
    return ggiPolyPatch_.master();
}


bool Foam::ggiFvPatch::fineLevel() const
{
    return true;
}


Foam::label Foam::ggiFvPatch::shadowIndex() const
{
    return ggiPolyPatch_.shadowIndex();
}


const Foam::ggiLduInterface& Foam::ggiFvPatch::shadowInterface() const
{
    const fvPatch& p = this->boundaryMesh()[ggiPolyPatch_.shadowIndex()];

    return refCast<const ggiLduInterface>(p);
}


Foam::label Foam::ggiFvPatch::zoneSize() const
{
    return ggiPolyPatch_.zone().size();
}


const Foam::labelList& Foam::ggiFvPatch::zoneAddressing() const
{
    return ggiPolyPatch_.zoneAddressing();
}


const Foam::labelListList& Foam::ggiFvPatch::addressing() const
{
    if (ggiPolyPatch_.master())
    {
        return ggiPolyPatch_.patchToPatch().masterAddr();
    }
    else
    {
        return ggiPolyPatch_.patchToPatch().slaveAddr();
    }
}


bool Foam::ggiFvPatch::localParallel() const
{
    return ggiPolyPatch_.localParallel();
}


const Foam::scalarListList& Foam::ggiFvPatch::weights() const
{
    if (ggiPolyPatch_.master())
    {
        return ggiPolyPatch_.patchToPatch().masterWeights();
    }
    else
    {
        return ggiPolyPatch_.patchToPatch().slaveWeights();
    }
}


Foam::tmp<Foam::labelField> Foam::ggiFvPatch::interfaceInternalField
(
    const unallocLabelList& internalData
) const
{
    return patchInternalField(internalData);
}


void Foam::ggiFvPatch::initTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& interfaceData
) const
{
    labelTransferBuffer_ = interfaceData;
}


Foam::tmp<Foam::labelField> Foam::ggiFvPatch::transfer
(
    const Pstream::commsTypes,
    const unallocLabelList& interfaceData
) const
{
    return this->shadow().labelTransferBuffer();
}


void Foam::ggiFvPatch::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{
    labelTransferBuffer_ = patchInternalField(iF);
}


Foam::tmp<Foam::labelField> Foam::ggiFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes,
    const unallocLabelList& iF
) const
{
    return shadow().labelTransferBuffer();
}



// ************************************************************************* //
