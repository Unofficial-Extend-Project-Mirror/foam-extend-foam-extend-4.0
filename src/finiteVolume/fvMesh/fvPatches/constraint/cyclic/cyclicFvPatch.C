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

#include "cyclicFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(cyclicFvPatch, 0);
addToRunTimeSelectionTable(fvPatch, cyclicFvPatch, polyPatch);


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Make patch weighting factors
void cyclicFvPatch::makeWeights(scalarField& w) const
{
    const scalarField& magFa = magSf();

    // Note: mag in the dot-product.
    // For all valid meshes, the non-orthogonality will be less that
    // 90 deg and the dot-product will be positive.  For invalid
    // meshes (d & s <= 0), this will stabilise the calculation
    // but the result will be poor.  HJ, 24/Aug/2011
    scalarField deltas = mag(nf() & fvPatch::delta());
    label sizeby2 = deltas.size()/2;

    scalar maxMatchError = 0;
    label errorFace = -1;

    for (label facei = 0; facei < sizeby2; facei++)
    {
        scalar avFa = (magFa[facei] + magFa[facei + sizeby2])/2.0;

        if
        (
            mag(magFa[facei] - magFa[facei + sizeby2])/avFa
          > polyPatch::matchTol_
        )
        {
            // Found error.  Look for largest matching error
            maxMatchError =
                Foam::max
                (
                    maxMatchError,
                    mag(magFa[facei] - magFa[facei + sizeby2])/avFa
                );

            errorFace = facei;
        }

        scalar di = deltas[facei];
        scalar dni = deltas[facei + sizeby2];

        w[facei] = dni/(di + dni);
        w[facei + sizeby2] = 1 - w[facei];
    }

    // Check for error in matching
    if (maxMatchError > polyPatch::matchTol_)
    {
        scalar avFa = (magFa[errorFace] + magFa[errorFace + sizeby2])/2.0;

        FatalErrorIn("cyclicFvPatch::makeWeights(scalarField& w) const")
            << "face " << errorFace << " and " << errorFace + sizeby2
            <<  " areas do not match by "
            << 100*mag(magFa[errorFace] - magFa[errorFace + sizeby2])/avFa
            << "% -- possible face ordering problem." << nl
            << "Cyclic area match tolerance = "
            << polyPatch::matchTol_ << " patch: " << name()
            << abort(FatalError);
    }
}


// Make patch face - neighbour cell distances
void cyclicFvPatch::makeDeltaCoeffs(scalarField& dc) const
{
    vectorField d = delta();
    vectorField n = nf();
    label sizeby2 = d.size()/2;

    for (label facei = 0; facei < sizeby2; facei++)
    {
        // Stabilised form for bad meshes.  HJ, 24/Aug/2011
        dc[facei] = 1.0/max(n[facei] & d[facei], 0.05*mag(d[facei]));
        dc[facei + sizeby2] = dc[facei];
    }
}


// Return delta (P to N) vectors across coupled patch
tmp<vectorField> cyclicFvPatch::delta() const
{
    vectorField patchD = fvPatch::delta();
    label sizeby2 = patchD.size()/2;

    tmp<vectorField> tpdv(new vectorField(patchD.size()));
    vectorField& pdv = tpdv();

    // To the transformation if necessary
    if (parallel())
    {
        for (label facei = 0; facei < sizeby2; facei++)
        {
            vector ddi = patchD[facei];
            vector dni = patchD[facei + sizeby2];

            pdv[facei] = ddi - dni;
            pdv[facei + sizeby2] = -pdv[facei];
        }
    }
    else
    {
        for (label facei = 0; facei < sizeby2; facei++)
        {
            vector ddi = patchD[facei];
            vector dni = patchD[facei + sizeby2];

            pdv[facei] = ddi - transform(forwardT()[0], dni);
            pdv[facei + sizeby2] = -transform(reverseT()[0], pdv[facei]);
        }
    }

    return tpdv;
}


tmp<labelField> cyclicFvPatch::interfaceInternalField
(
    const unallocLabelList& internalData
) const
{
    return patchInternalField(internalData);
}


tmp<labelField> cyclicFvPatch::transfer
(
    const Pstream::commsTypes,
    const unallocLabelList& interfaceData
) const
{
    tmp<labelField> tpnf(new labelField(this->size()));
    labelField& pnf = tpnf();

    label sizeby2 = this->size()/2;

    for (label facei=0; facei<sizeby2; facei++)
    {
        pnf[facei] = interfaceData[facei + sizeby2];
        pnf[facei + sizeby2] = interfaceData[facei];
    }

    return tpnf;
}


tmp<labelField> cyclicFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{
    const unallocLabelList& faceCells = this->patch().faceCells();

    tmp<labelField> tpnf(new labelField(this->size()));
    labelField& pnf = tpnf();

    label sizeby2 = this->size()/2;

    for (label facei=0; facei<sizeby2; facei++)
    {
        pnf[facei] = iF[faceCells[facei + sizeby2]];
        pnf[facei + sizeby2] = iF[faceCells[facei]];
    }

    return tpnf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
