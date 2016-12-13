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
    Harmonic-mean differencing scheme class.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceScalarField> Foam::harmonic<Type>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& phi
) const
{
    // Implement weights-based stabilised harmonic interpolation using
    // magnitude of type
    // Algorithm:
    // 1) calculate magnitude of internal field and neighbour field
    // 2) calculate harmonic mean magnitude
    // 3) express harmonic mean magnitude as: mean = w*mOwn + (1 - w)*mNei
    // 4) Based on above, calculate w = (mean - mNei)/(mOwn - mNei)
    // 5) Use weights to interpolate values

    tmp<surfaceScalarField> tw
    (
        new surfaceScalarField
        (
            IOobject
            (
                "harmonicWeightingFactors" + phi.name(),
                this->mesh().time().timeName(),
                this->mesh()
            ),
            this->mesh() ,
            dimless
        )
    );

    const surfaceScalarField& deltaCoeffs = this->mesh().deltaCoeffs();
    const surfaceScalarField& weights = this->mesh().weights();

    const magLongDelta& mld = magLongDelta::New(this->mesh());
    const surfaceScalarField& longDelta = mld.magDelta();

    surfaceScalarField& w = tw();

    const unallocLabelList& owner = this->mesh().owner();
    const unallocLabelList& neighbour = this->mesh().neighbour();

    scalarField magPhi = mag(phi);

    scalarField& wIn = w.internalField();

    // Larger small for complex arithmetic accuracy
    const scalar kSmall = 1000*SMALL;

    // Calculate internal weights using field magnitude
    forAll (owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        scalar mOwn = magPhi[own]/(1 - weights[faceI]);
        scalar mNei = magPhi[nei]/weights[faceI];

        scalar den = magPhi[nei] - magPhi[own];

        scalar mean = mOwn*mNei/
            ((mOwn + mNei)*longDelta[faceI]*deltaCoeffs[faceI]);

        // Note: complex arithmetic requires extra accuracy
        // This is a division of two close subtractions
        // HJ, 28/Sep/2011
        if (mag(den) > kSmall)
        {
            // Limit weights for round-off safety
            wIn[faceI] =
                Foam::max(0, Foam::min((magPhi[nei] - mean)/den, 1));
        }
        else
        {
            wIn[faceI] = 0.5;
        }
    }

    forAll (phi.boundaryField(), pi)
    {
        fvsPatchScalarField& wp = w.boundaryField()[pi];

        const fvPatchField<Type>& pPhi = phi.boundaryField()[pi];

        if (pPhi.coupled())
        {
            wp = this->weights
            (
                pPhi.patchInternalField(),
                pPhi.patchNeighbourField(),
                pPhi.patch()
            );
        }
        else
        {
            // Boundary weights for uncoupled patches are 1
            wp = 1;
        }
    }

    return tw;
}


template<class Type>
Foam::tmp<Foam::scalarField> Foam::harmonic<Type>::weights
(
    const Field<Type>& fOwn,
    const Field<Type>& fNei,
    const fvPatch& patch
) const
{
    // Implement weights-based stabilised harmonic interpolation using
    // magnitude of type
    // Algorithm:
    // 1) calculate magnitude of internal field and neighbour field
    // 2) calculate harmonic mean magnitude
    // 3) express harmonic mean magnitude as: mean = w*mOwn + (1 - w)*mNei
    // 4) Based on above, calculate w = (mean - mNei)/(mOwn - mNei)
    // 5) Use weights to interpolate values

    tmp<scalarField> tweights(new scalarField(fOwn.size(), 0.5));
    scalarField& weights = tweights();

    // Larger small for complex arithmetic accuracy
    const scalar kSmall = 1000*SMALL;

    // Mag long deltas are identical on both sides.  HJ, 28/Sep/2011
    const magLongDelta& mld = magLongDelta::New(this->mesh());

    scalarField magPhiOwn = mag(fOwn);
    scalarField magPhiNei = mag(fNei);

    const scalarField& pWeights = patch.weights();
    const scalarField& pDeltaCoeffs = patch.deltaCoeffs();
    const scalarField& pLongDelta = mld.magDelta(patch.index());

    forAll (weights, faceI)
    {
        scalar mOwn = magPhiOwn[faceI]/(1 - pWeights[faceI]);
        scalar mNei = magPhiNei[faceI]/pWeights[faceI];

        scalar den = magPhiNei[faceI] - magPhiOwn[faceI];

        // Note: complex arithmetic requires extra accuracy
        // This is a division of two close subtractions
        // HJ, 28/Sep/2011
        if (mag(den) > kSmall)
        {
            scalar mean = mOwn*mNei/
                (
                    (mOwn + mNei)*
                    pLongDelta[faceI]*
                    pDeltaCoeffs[faceI]
                );

            // Limit weights for round-off safety
            weights[faceI] =
                Foam::max(0, Foam::min((magPhiNei[faceI] - mean)/den, 1));
        }
        else
        {
            // Use 0.5 weights
        }
    }

    return tweights;
}


// ************************************************************************* //
