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
