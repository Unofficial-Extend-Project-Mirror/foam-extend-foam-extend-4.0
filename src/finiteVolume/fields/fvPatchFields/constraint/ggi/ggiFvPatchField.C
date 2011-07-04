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

// Return neighbour field
template<class Type>
tmp<Field<Type> > ggiFvPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->internalField();

    // Get shadow face-cells and assemble shadow field
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
        vectorField nHat = this->patch().nf();

        // Use mirrored neighbour field for interpolation
        // HJ, 21/Jan/2009
        Field<Type> bridgeField =
            transform(I - 2.0*sqr(nHat), this->patchInternalField());

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

    if (ggiPatch_.bridgeOverlap())
    {
        // Symmetry treatment used for overlap
        vectorField nHat = this->patch().nf();

        Field<Type> bridgeField =
        (
            this->patchInternalField()
          + transform(I - 2.0*sqr(nHat), this->patchInternalField())
        )/2.0;

        ggiPatch_.bridge(bridgeField, pf);
    }

    Field<Type>::operator=(pf);
}


template<class Type>
void ggiFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
}


template<class Type>
void ggiFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    // Communication is allowed either before or after processor
    // patch comms.  HJ, 11/Jul/2011

    // Get shadow face-cells and assemble shadow field
    const unallocLabelList& sfc = ggiPatch_.shadow().faceCells();

    scalarField sField(sfc.size());

    forAll (sField, i)
    {
        sField[i] = psiInternal[sfc[i]];
    }

    scalarField pnf = ggiPatch_.interpolate(sField);

    // Multiply the field by coefficients and add into the result
    const unallocLabelList& fc = ggiPatch_.faceCells();

    forAll(fc, elemI)
    {
        result[fc[elemI]] -= coeffs[elemI]*pnf[elemI];
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
    const Pstream::commsTypes
) const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
