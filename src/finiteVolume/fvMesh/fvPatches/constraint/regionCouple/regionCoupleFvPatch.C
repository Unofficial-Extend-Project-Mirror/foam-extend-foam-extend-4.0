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
    Region couple patch

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "regionCoupleFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "fvBoundaryMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionCoupleFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, regionCoupleFvPatch, polyPatch);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Make patch weighting factors
void Foam::regionCoupleFvPatch::makeWeights(scalarField& w) const
{
    if (rcPolyPatch_.attached())
    {
        vectorField n = nf();
        scalarField nfc = n & (rcPolyPatch_.reconFaceCellCentres() - Cf());

        w = nfc/((n & (Cf() - Cn())) + nfc);
    }
    else
    {
        fvPatch::makeWeights(w);
    }
}


// Make patch face - neighbour cell distances
void Foam::regionCoupleFvPatch::makeDeltaCoeffs(scalarField& dc) const
{
    if (rcPolyPatch_.attached())
    {
        dc = (1.0 - weights())/(nf() & fvPatch::delta());
    }
    else
    {
        fvPatch::makeDeltaCoeffs(dc);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

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


// Return delta (P to N) vectors across coupled patch
Foam::tmp<Foam::vectorField> Foam::regionCoupleFvPatch::delta() const
{
    if (rcPolyPatch_.attached())
    {
        return rcPolyPatch_.reconFaceCellCentres() - Cn();
    }
    else
    {
        return fvPatch::delta();
    }
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
    transferBuffer_ = interfaceData;
}


Foam::tmp<Foam::labelField> Foam::regionCoupleFvPatch::transfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& interfaceData
) const
{
    //HJ  Should this be mapped?  22/Jun/2007
    return shadow().transferBuffer();
}


void Foam::regionCoupleFvPatch::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{
    transferBuffer_ = patchInternalField(iF);
}


Foam::tmp<Foam::labelField> Foam::regionCoupleFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{
    //HJ  Should this be mapped?  22/Jun/2007
    return shadow().transferBuffer();
}


// ************************************************************************* //
