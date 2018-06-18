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
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Contributor
    Martin Beaudoin, Hydro-Quebec, (2008)

Note on parallelisation
    In order to handle parallelisation correctly, I need to rely on the fact
    that all patches that require a global gather-scatter come before
    processor patches.  In that case, the communication pattern
    will be correct without intervention.  HJ, 6/Aug/2009

\*---------------------------------------------------------------------------*/

#include "ggiFvPatchField.H"
#include "symmTransformField.H"
#include "coeffFields.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
ggiFvPatchField<Type>::ggiFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    ggiPatch_(refCast<const ggiFvPatch>(p))
{}


template<class Type>
ggiFvPatchField<Type>::ggiFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict, false),
    ggiPatch_(refCast<const ggiFvPatch>(p))
{
    if (!isType<ggiFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "ggiFvPatchField<Type>::ggiFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not ggi type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }

    if (!dict.found("value"))
    {
        // Grab the internal value for initialisation. (?) HJ, 27/Feb/2009
        fvPatchField<Type>::operator=(this->patchInternalField()());
    }
}


template<class Type>
ggiFvPatchField<Type>::ggiFvPatchField
(
    const ggiFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    ggiPatch_(refCast<const ggiFvPatch>(p))
{
    if (!isType<ggiFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "ggiFvPatchField<Type>::ggiFvPatchField\n"
            "(\n"
            "    const ggiFvPatchField<Type>& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
ggiFvPatchField<Type>::ggiFvPatchField
(
    const ggiFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    ggiLduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    ggiPatch_(refCast<const ggiFvPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const ggiFvPatchField<Type>& ggiFvPatchField<Type>::shadowPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
    static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
    (
        this->internalField()
    );

    return refCast<const ggiFvPatchField<Type> >
    (
        fld.boundaryField()[ggiPatch_.shadowIndex()]
    );
}


template<class Type>
tmp<Field<Type> > ggiFvPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->internalField();

    // Get shadow face-cells and assemble shadow field
    // This is a patchInternalField of neighbour but access is inconvenient.
    // Assemble by hand. HJ, 27/Sep/2011
    const unallocLabelList& sfc = ggiPatch_.shadow().faceCells();

    Field<Type> sField(sfc.size());

    forAll (sField, i)
    {
        sField[i] = iField[sfc[i]];
    }

    tmp<Field<Type> > tpnf(ggiPatch_.interpolate(sField));
    Field<Type>& pnf = tpnf();

    if (ggiPatch_.bridgeOverlap())
    {
        // Use mirrored neighbour field for interpolation. Note: mirroring needs
        // to take into account the weights, i.e. how "far" we are actually
        // mirroring. VV, 19/Jan/2018.
        const Field<Type> mirrorField =
        transform
        (
            (I - sqr(this->patch().nf())/(1.0 - ggiPatch_.fvPatch::weights())),
            this->patchInternalField()
        );

        // Set mirror values to fully uncovered faces
        ggiPatch_.setUncoveredFaces(mirrorField, pnf);

        // Add part of the mirror field to partially covered faces
        ggiPatch_.addToPartialFaces(mirrorField, pnf);
    }

    return tpnf;
}


template<class Type>
void ggiFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type> pf
    (
        this->patch().weights()*this->patchInternalField()
      + (1.0 - this->patch().weights())*this->patchNeighbourField()
    );

    // Note: bridging and correction of partially overlapping faces taken into
    // account in patchNeighbourField(). VV, 16/Oct/2017.

    Field<Type>::operator=(pf);
}


template<class Type>
void ggiFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsTypes
)
{
    coupledFvPatchField<Type>::evaluate(commsTypes);
}


template<class Type>
void ggiFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    // Communication is allowed either before or after processor
    // patch comms.  HJ, 11/Jul/2011

    // Get shadow face-cells and assemble shadow field
    // Consider using shadowField().patchInternalField() ? HJ, 18/Feb/2016
    const unallocLabelList& sfc = ggiPatch_.shadow().faceCells();

    scalarField sField(sfc.size());

    forAll (sField, i)
    {
        sField[i] = psiInternal[sfc[i]];
    }

    scalarField pnf = ggiPatch_.interpolate(sField);

    if (ggiPatch_.bridgeOverlap())
    {
        // Note: this implicit treatment does not really work implicitly for
        // types with rank > 0 (everything above scalar) if the symmetry plane
        // is not aligned with one of the coordinate axes. VV, 18/Oct/2017.

        // Use mirrored neighbour field for interpolation. Note: mirroring needs
        // to take into account the weights, i.e. how "far" we are actually
        // mirroring. VV, 19/Jan/2018.
        const scalarField mirrorField =
        transform
        (
            (I - sqr(this->patch().nf())/(1.0 - ggiPatch_.fvPatch::weights())),
            ggiPatch_.patchInternalField(psiInternal)
        );

        // Only set fully uncovered faces. Partially covered faces taken into
        // account by manipulating value and gradient matrix coefficients
        ggiPatch_.setUncoveredFaces(mirrorField, pnf);
    }

    // Multiply the field by coefficients and add into the result
    const unallocLabelList& fc = ggiPatch_.faceCells();

    if (switchToLhs)
    {
        forAll(fc, elemI)
        {
            result[fc[elemI]] += coeffs[elemI]*pnf[elemI];
        }
    }
    else
    {
        forAll(fc, elemI)
        {
            result[fc[elemI]] -= coeffs[elemI]*pnf[elemI];
        }
    }
}


template<class Type>
void ggiFvPatchField<Type>::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes,
    const bool switchToLhs
) const
{}


template<class Type>
void ggiFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const Field<Type>& psiInternal,
    Field<Type>& result,
    const BlockLduMatrix<Type>&,
    const CoeffField<Type>& coeffs,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    // Communication is allowed either before or after processor
    // patch comms.  HJ, 11/Jul/2011

    // Get shadow face-cells and assemble shadow patch internal field
    // Consider using shadowField().patchInternalField() ? HJ, 18/Feb/2016
    const unallocLabelList& sfc = ggiPatch_.shadow().faceCells();

    Field<Type> sField(sfc.size());

    forAll (sField, i)
    {
        sField[i] = psiInternal[sfc[i]];
    }

    Field<Type> pnf = ggiPatch_.interpolate(sField);

    if (ggiPatch_.bridgeOverlap())
    {
        // Use mirrored neighbour field for interpolation. Note: mirroring needs
        // to take into account the weights, i.e. how "far" we are actually
        // mirroring. VV, 19/Jan/2018.
        const Field<Type> mirrorField =
        transform
        (
            (I - sqr(this->patch().nf())/(1.0 - ggiPatch_.fvPatch::weights())),
            ggiPatch_.patchInternalField(psiInternal)
        );

        // Only set fully uncovered faces. Partially covered faces taken into
        // account by manipulating value and gradient matrix coefficients
        ggiPatch_.setUncoveredFaces(mirrorField, pnf);
    }

    // Multiply neighbour field with coeffs and re-use pnf for result
    // of multiplication
    multiply(pnf, coeffs, pnf);

    const unallocLabelList& fc = ggiPatch_.faceCells();

    if (switchToLhs)
    {
        forAll (fc, elemI)
        {
            result[fc[elemI]] += pnf[elemI];
        }
    }
    else
    {
        forAll (fc, elemI)
        {
            result[fc[elemI]] -= pnf[elemI];
        }
    }
}


template<class Type>
void ggiFvPatchField<Type>::updateInterfaceMatrix
(
    const Field<Type>&,
    Field<Type>&,
    const BlockLduMatrix<Type>&,
    const CoeffField<Type>&,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    // Does nothing
}


template<class Type>
void Foam::ggiFvPatchField<Type>::patchFlux
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& flux,
    const fvMatrix<Type>& matrix
) const
{
    // Since we have adjusted the internal/boundary coefficients in the
    // manipulateMatrix member function below, we must not use
    // patchNeighbourField for reconstructing the flux. We only need to use
    // interpolated shadow field. VV, 6/Mar/2018.

    // Get patch ID
    const label patchI = this->patch().index();

    // Get internal field
    const Field<Type>& iField = this->internalField();

    // Get shadow face-cells and assemble shadow field
    const unallocLabelList& sfc = ggiPatch_.shadow().faceCells();

    Field<Type> sField(sfc.size());
    forAll (sField, i)
    {
        sField[i] = iField[sfc[i]];
    }

    // Interpolate shadow to this side. Note: must not bridge since internal
    // coeffs and boundary coeffs take it into account
    Field<Type> neighbourField(ggiPatch_.interpolate(sField));

    if (ggiPatch_.bridgeOverlap())
    {
        const Field<Type> mirrorField =
        transform
        (
            (I - sqr(this->patch().nf())/(1.0 - ggiPatch_.fvPatch::weights())),
            this->patchInternalField()
        );

        // Only set fully uncovered faces. Partially covered faces taken into
        // account by manipulating value and gradient matrix coefficients
        ggiPatch_.setUncoveredFaces(mirrorField, neighbourField);
    }

    // Calculate the flux with correct neighbour field (fully uncovered faces
    // bridged, while partially uncovered faces taken into account by
    // manipulating value and gradient matrix coefficients in order to ensure
    // conservation for both convection and diffusion part across partially
    // overlapping faces). VV, 14/Mar/2018.
    flux.boundaryField()[patchI] =
        cmptMultiply
        (
            matrix.internalCoeffs()[patchI],
            this->patchInternalField()
        )
      - cmptMultiply
        (
            matrix.boundaryCoeffs()[patchI],
            neighbourField
        );
}


template<class Type>
void Foam::ggiFvPatchField<Type>::manipulateValueCoeffs
(
    fvMatrix<Type>& matrix
) const
{
    // Conservative treatment for convection across bridged overlap
    // Note: master corrects both sets of coefficients (master and slave)
    if (ggiPatch_.bridgeOverlap() && ggiPatch_.master())
    {
        // Get this patch and shadow patch index. In this case, this patch is
        // always master and shadow patch is always slave (second condition in
        // the if statement)
        const label mPatchI = this->patch().index();
        const label sPatchI = ggiPatch_.shadowIndex();

        // Get all matrix coefficients
        Field<Type>& masterIC = matrix.internalCoeffs()[mPatchI];
        Field<Type>& masterBC = matrix.boundaryCoeffs()[mPatchI];
        Field<Type>& slaveIC = matrix.internalCoeffs()[sPatchI];
        Field<Type>& slaveBC = matrix.boundaryCoeffs()[sPatchI];

        // Get surface area magnitudes
        const scalarField& magSfMaster = ggiPatch_.magSf();
        const scalarField& magSfSlave = ggiPatch_.shadow().magSf();

        // Interpolate all coeffs from each side on the other side. Note: need
        // to normalize by surface areas since it is "hidden" inside the
        // convection coefficient (usually Sf & Uf), before interpolation. After
        // interpolation, need to multiply with surface areas on this side.

        // Interpolation from slave to this (master) side
        const Field<Type> masterInterpolatedIC =
            magSfSlave*ggiPatch_.shadow().interpolate(masterIC/magSfMaster);
        const Field<Type> masterInterpolatedBC =
            magSfSlave*ggiPatch_.shadow().interpolate(masterBC/magSfMaster);

        // Interpolation from this (master) side to slave
        const Field<Type> slaveInterpolatedIC =
            magSfMaster*ggiPatch_.interpolate(slaveIC/magSfSlave);
        const Field<Type> slaveInterpolatedBC =
            magSfMaster*ggiPatch_.interpolate(slaveBC/magSfSlave);


        // Set partially covered master coeffs using slave data
        ggiPatch_.setPartialFaces(slaveInterpolatedBC, masterIC);
        ggiPatch_.setPartialFaces(slaveInterpolatedIC, masterBC);

        // Scale partially overlapping boundary coeffs to ensure conservation
        ggiPatch_.scalePartialFaces(masterBC);


        // Set partially covered slave coeffs using master data
        ggiPatch_.shadow().setPartialFaces(masterInterpolatedBC, slaveIC);
        ggiPatch_.shadow().setPartialFaces(masterInterpolatedIC, slaveBC);

        // Scale partially overlapping boundary coeffs to ensure conservation
        ggiPatch_.shadow().scalePartialFaces(slaveBC);
    }
}


template<class Type>
void Foam::ggiFvPatchField<Type>::manipulateGradientCoeffs
(
    fvMatrix<Type>& matrix
) const
{
    // Conservative treatment for convection across bridged overlap
    // Note: master corrects both sets of coefficients (master and slave)
    if (ggiPatch_.bridgeOverlap() && ggiPatch_.master())
    {
        // Get this patch and shadow patch index. In this case, this patch is
        // always master and shadow patch is always slave (second condition in
        // the if statement)
        const label mPatchI = this->patch().index();
        const label sPatchI = ggiPatch_.shadowIndex();

        // Get all matrix coefficients
        Field<Type>& masterIC = matrix.internalCoeffs()[mPatchI];
        Field<Type>& masterBC = matrix.boundaryCoeffs()[mPatchI];
        Field<Type>& slaveIC = matrix.internalCoeffs()[sPatchI];
        Field<Type>& slaveBC = matrix.boundaryCoeffs()[sPatchI];

        // Get surface area magnitudes
        const scalarField& magSfMaster = ggiPatch_.magSf();
        const scalarField& magSfSlave = ggiPatch_.shadow().magSf();

        // Interpolate all coeffs from each side on the other side. Note: need
        // to normalize by surface areas since it is "hidden" inside the
        // diffusion coefficient (usually |Sf|/|df|), before interpolation.
        // After interpolation, need to multiply with surface areas on this
        // side.

        // Interpolation from slave to this (master) side
        const Field<Type> masterInterpolatedIC =
            magSfSlave*ggiPatch_.shadow().interpolate(masterIC/magSfMaster);
        const Field<Type> masterInterpolatedBC =
            magSfSlave*ggiPatch_.shadow().interpolate(masterBC/magSfMaster);

        const Field<Type> slaveInterpolatedIC =
            magSfMaster*ggiPatch_.interpolate(slaveIC/magSfSlave);
        const Field<Type> slaveInterpolatedBC =
            magSfMaster*ggiPatch_.interpolate(slaveBC/magSfSlave);


        // Set partially covered master coeffs using slave data
        ggiPatch_.setPartialFaces(slaveInterpolatedBC, masterIC);
        ggiPatch_.setPartialFaces(slaveInterpolatedIC, masterBC);

        // Scale partially overlapping master coeffs to ensure conservation
        ggiPatch_.scalePartialFaces(masterIC);
        ggiPatch_.scalePartialFaces(masterBC);

        // Add to partially overlapping master internal coeffs to ensure
        // conservation. Note: boundary coeffs must be negated
        ggiPatch_.addToPartialFaces((-masterBC)(), masterIC);


        // Set partially covered slave coeffs using master data
        ggiPatch_.shadow().setPartialFaces(masterInterpolatedBC, slaveIC);
        ggiPatch_.shadow().setPartialFaces(masterInterpolatedIC, slaveBC);

        // Scale partially overlapping slave coeffs to ensure conservation
        ggiPatch_.shadow().scalePartialFaces(slaveIC);
        ggiPatch_.shadow().scalePartialFaces(slaveBC);

        // Add to partially overlapping slave internal coeffs to ensure
        // conservation. Note: boundary coeffs must be negated
        ggiPatch_.shadow().addToPartialFaces((-slaveBC)(), slaveIC);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
