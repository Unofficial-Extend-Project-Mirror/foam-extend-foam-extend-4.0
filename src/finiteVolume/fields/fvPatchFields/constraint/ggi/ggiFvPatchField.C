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
        // Symmetry treatment used for overlap
        const vectorField nHat = this->patch().nf();

        // Use mirrored neighbour field for interpolation
        // HJ, 21/Jan/2009
        const Field<Type> bridgeField =
            transform(I - 2.0*sqr(nHat), this->patchInternalField());

        // Note: bridging now takes into account fully uncovered and partially
        // covered faces. VV, 18/Oct/2017.
        ggiPatch_.bridge(bridgeField, pnf);
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
    const Pstream::commsTypes
)
{
    fvPatchField<Type>::evaluate();
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
        // Note: will not work properly for types with rank > 0 (everything
        // above scalar) if the symmetry plane is not aligned with one of the
        // coordinate axes. VV, 18/Oct/2017.
        const scalarField bridgeField =
            transform
            (
                I - 2.0*sqr(this->patch().nf()),
                ggiPatch_.patchInternalField(psiInternal)
            );

        ggiPatch_.bridge(bridgeField, pnf);
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
        // Note: will not work properly for types with rank > 0 (everything
        // above scalar) if the symmetry plane is not aligned with one of the
        // coordinate axes. VV, 18/Oct/2017.
        const Field<Type> bridgeField =
            transform
            (
                I - 2.0*sqr(this->patch().nf()),
                ggiPatch_.patchInternalField(psiInternal)
            );

        ggiPatch_.bridge(bridgeField, pnf);
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
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
