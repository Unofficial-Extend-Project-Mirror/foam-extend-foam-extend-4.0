/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

Description
    Region couple patch

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "regionCoupleFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "fvBoundaryMesh.H"
#include "foamTime.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionCoupleFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, regionCoupleFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionCoupleFvPatch::~regionCoupleFvPatch()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Make patch weighting factors
void Foam::regionCoupleFvPatch::makeWeights(scalarField& w) const
{
    if (rcPolyPatch_.coupled())
    {
        if (rcPolyPatch_.master())
        {
            vectorField n = nf();

            // Note: mag in the dot-product.
            // For all valid meshes, the non-orthogonality will be less than
            // 90 deg and the dot-product will be positive.  For invalid
            // meshes (d & s <= 0), this will stabilise the calculation
            // but the result will be poor.  HJ, 24/Aug/2011
            scalarField nfc =
                mag(n & (rcPolyPatch_.reconFaceCellCentres() - Cf()));

            w = nfc/(mag(n & (Cf() - Cn())) + nfc);

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

            scalarField oneMinusW = 1 - masterWeights;

            w = interpolate(oneMinusW);

            if (bridgeOverlap())
            {
                // Set overlap weights to 0.5 and use mirrored neighbour field
                // for interpolation.  HJ, 21/Jan/2009
                bridge(scalarField(size(), 0.5), w);
            }
        }
    }
    else
    {
        fvPatch::makeWeights(w);
    }
}


// Make patch face - neighbour cell distances
void Foam::regionCoupleFvPatch::makeDeltaCoeffs(scalarField& dc) const
{
    if (rcPolyPatch_.coupled())
    {
        if (rcPolyPatch_.master())
        {
            // Stabilised form for bad meshes.  HJ, 24/Aug/2011
            vectorField d = delta();

            dc = 1.0/max(nf() & d, 0.05*mag(d));

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
    else
    {
        fvPatch::makeDeltaCoeffs(dc);
    }
}


// Make patch face non-orthogonality correction vectors
void Foam::regionCoupleFvPatch::makeCorrVecs(vectorField& cv) const
{
    if (rcPolyPatch_.coupled())
    {
        // Non-orthogonality correction in attached state identical to ggi
        // interface

        // Calculate correction vectors on coupled patches
        const scalarField& patchDeltaCoeffs = deltaCoeffs();

        vectorField patchDeltas = delta();
        vectorField n = nf();

        // If non-orthogonality is over 90 deg, kill correction vector
        // HJ, 6/Jan/2011
        cv = pos(patchDeltas & n)*(n - patchDeltas*patchDeltaCoeffs);
    }
    else
    {
        // No correction in detached state.  HJ, 26/Jul/2011
        cv = vector::zero;
    }
}


// Return delta (P to N) vectors across coupled patch
Foam::tmp<Foam::vectorField> Foam::regionCoupleFvPatch::delta() const
{
    if (rcPolyPatch_.coupled())
    {
        if (rcPolyPatch_.master())
        {
            tmp<vectorField> tDelta =
                rcPolyPatch_.reconFaceCellCentres() - Cn();

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
                shadow().Cn() - rcPolyPatch_.shadow().reconFaceCellCentres()
            );

            if (bridgeOverlap())
            {
                vectorField bridgeDeltas = Cf() - Cn();

                bridge(bridgeDeltas, tDelta());
            }

            return tDelta;
        }
    }
    else
    {
        return fvPatch::delta();
    }
}


bool Foam::regionCoupleFvPatch::coupled() const
{
    return rcPolyPatch_.coupled();
}


const Foam::fvMesh& Foam::regionCoupleFvPatch::shadowRegion() const
{
    return
        boundaryMesh().mesh().objectRegistry::parent().lookupObject<fvMesh>
        (
            rcPolyPatch_.shadowRegionName()
        );
}


const Foam::regionCoupleFvPatch& Foam::regionCoupleFvPatch::shadow() const
{
    const fvPatch& p =
        shadowRegion().boundary()[rcPolyPatch_.shadowIndex()];

    return refCast<const regionCoupleFvPatch>(p);
}


bool Foam::regionCoupleFvPatch::master() const
{
    return rcPolyPatch_.master();
}


bool Foam::regionCoupleFvPatch::fineLevel() const
{
    return true;
}


Foam::label Foam::regionCoupleFvPatch::shadowIndex() const
{
    return rcPolyPatch_.shadowIndex();
}


const Foam::ggiLduInterface&
Foam::regionCoupleFvPatch::shadowInterface() const
{
    const fvPatch& p =
        shadowRegion().boundary()[rcPolyPatch_.shadowIndex()];

    return refCast<const ggiLduInterface>(p);
}


Foam::label Foam::regionCoupleFvPatch::interfaceSize() const
{
    return rcPolyPatch_.size();
}


Foam::label Foam::regionCoupleFvPatch::zoneSize() const
{
    return rcPolyPatch_.zone().size();
}


const Foam::labelList& Foam::regionCoupleFvPatch::zoneAddressing() const
{
    return rcPolyPatch_.zoneAddressing();
}


const Foam::labelListList& Foam::regionCoupleFvPatch::addressing() const
{
    if (rcPolyPatch_.master())
    {
        return rcPolyPatch_.patchToPatch().masterAddr();
    }
    else
    {
        return rcPolyPatch_.patchToPatch().slaveAddr();
    }
}


bool Foam::regionCoupleFvPatch::localParallel() const
{
    return rcPolyPatch_.localParallel();
}


const Foam::mapDistribute& Foam::regionCoupleFvPatch::map() const
{
    return rcPolyPatch_.map();
}


const Foam::scalarListList& Foam::regionCoupleFvPatch::weights() const
{
    if (rcPolyPatch_.master())
    {
        return rcPolyPatch_.patchToPatch().masterWeights();
    }
    else
    {
        return rcPolyPatch_.patchToPatch().slaveWeights();
    }
}


void Foam::regionCoupleFvPatch::expandAddrToZone(labelField& lf) const
{
    // Missing code.  Activate for AMG solvers across regionCoupleFvPatch
    notImplemented
    (
        "void regionCoupleFvPatch::expandAddrToZone(labelField& lf) const"
    );
}


Foam::tmp<Foam::labelField> Foam::regionCoupleFvPatch::interfaceInternalField
(
    const unallocLabelList& internalData
) const
{
    return patchInternalField(internalData);
}


void Foam::regionCoupleFvPatch::initTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& interfaceData
) const
{
    labelTransferBuffer_ = interfaceData;
}


Foam::tmp<Foam::labelField> Foam::regionCoupleFvPatch::transfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& interfaceData
) const
{

    return shadow().labelTransferBuffer();
}


void Foam::regionCoupleFvPatch::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{
    labelTransferBuffer_ = patchInternalField(iF);
}


Foam::tmp<Foam::labelField> Foam::regionCoupleFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{
    return shadow().labelTransferBuffer();
}


// ************************************************************************* //
