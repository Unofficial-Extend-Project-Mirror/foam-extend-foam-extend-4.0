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
    Martin Beaudoin, Hydro-Quebec, (2008)

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "cyclicGgiFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicGgiFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, cyclicGgiFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Make patch weighting factors
void Foam::cyclicGgiFvPatch::makeWeights(scalarField& w) const
{
    // Calculation of weighting factors is performed from the master
    // position, using reconstructed shadow cell centres
    // HJ, 2/Aug/2007
    if (cyclicGgiPolyPatch_.master())
    {
        vectorField n = nf();
        scalarField nfc =
            n & (cyclicGgiPolyPatch_.reconFaceCellCentres() - Cf());

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
void Foam::cyclicGgiFvPatch::makeDeltaCoeffs(scalarField& dc) const
{
    if (cyclicGgiPolyPatch_.master())
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
        tmp<vectorField> tDelta =
            cyclicGgiPolyPatch_.reconFaceCellCentres() - Cn();

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
            shadow().Cn() - cyclicGgiPolyPatch_.shadow().reconFaceCellCentres()
        );

        if (bridgeOverlap())
        {
            vectorField bridgeDeltas = Cf() - Cn();

            bridge(bridgeDeltas, tDelta());
        }

        return tDelta;
    }
}


// ************************************************************************* //
