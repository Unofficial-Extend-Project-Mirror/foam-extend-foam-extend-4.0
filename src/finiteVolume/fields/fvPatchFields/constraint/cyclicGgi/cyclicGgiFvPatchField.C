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

Contributor:
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "cyclicGgiFvPatchField.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
cyclicGgiFvPatchField<Type>::cyclicGgiFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    cyclicGgiPatch_(refCast<const cyclicGgiFvPatch>(p))
{}


template<class Type>
cyclicGgiFvPatchField<Type>::cyclicGgiFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict, false),
    cyclicGgiPatch_(refCast<const cyclicGgiFvPatch>(p))
{
    if (!isType<cyclicGgiFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "cyclicGgiFvPatchField<Type>::cyclicGgiFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not cyclicGgi type. "
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
cyclicGgiFvPatchField<Type>::cyclicGgiFvPatchField
(
    const cyclicGgiFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    cyclicGgiPatch_(refCast<const cyclicGgiFvPatch>(p))
{
    if (!isType<cyclicGgiFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "cyclicGgiFvPatchField<Type>::cyclicGgiFvPatchField\n"
            "(\n"
            "    const cyclicGgiFvPatchField<Type>& ptf,\n"
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
cyclicGgiFvPatchField<Type>::cyclicGgiFvPatchField
(
    const cyclicGgiFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    ggiLduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    cyclicGgiPatch_(refCast<const cyclicGgiFvPatch>(ptf.patch()))
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return neighbour field
template<class Type>
tmp<Field<Type> > cyclicGgiFvPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->internalField();

    // Get shadow face-cells and assemble shadow field
    const unallocLabelList& sfc = cyclicGgiPatch_.shadow().faceCells();

    Field<Type> sField(sfc.size());

    forAll (sField, i)
    {
        sField[i] = iField[sfc[i]];
    }

    // Transformation is handled in interpolation.  HJ, 7/Jan/2009
    tmp<Field<Type> > tpnf(cyclicGgiPatch_.interpolate(sField));
    Field<Type>& pnf = tpnf();

    if (cyclicGgiPatch_.bridgeOverlap())
    {
        // Symmetry treatment used for overlap
        vectorField nHat = this->patch().nf();

        // Use mirrored internal field for neighbour
        // HJ, 27/Jan/2009
        Field<Type> bridgeField =
            transform(I - 2.0*sqr(nHat), this->patchInternalField());

        cyclicGgiPatch_.bridge(bridgeField, pnf);
    }

    return tpnf;
}


template<class Type>
void cyclicGgiFvPatchField<Type>::initEvaluate
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

    if (cyclicGgiPatch_.bridgeOverlap())
    {
        // Symmetry treatment used for overlap
        vectorField nHat = this->patch().nf();

        Field<Type> bridgeField =
        (
            this->patchInternalField()
          + transform(I - 2.0*sqr(nHat), this->patchInternalField())
        )/2.0;

        cyclicGgiPatch_.bridge(bridgeField, pf);
    }

    Field<Type>::operator=(pf);
}


template<class Type>
void cyclicGgiFvPatchField<Type>::evaluate
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
void cyclicGgiFvPatchField<Type>::initInterfaceMatrixUpdate
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
    const unallocLabelList& sfc = cyclicGgiPatch_.shadow().faceCells();

    scalarField sField(sfc.size());

    forAll (sField, i)
    {
        sField[i] = psiInternal[sfc[i]];
    }

    // Transform according to the transformation tensor, using the slave
    // side transform.  Warning: forwardT() corresponds to the slave
    // patch.  HJ, 12/Jan/2009
    transformCoupleField(sField, cmpt);

    // Note: scalar interpolate does not get a transform, so this is safe
    // HJ, 12/Jan/2009
    scalarField pnf = cyclicGgiPatch_.interpolate(sField);

    // Multiply the field by coefficients and add into the result
    const unallocLabelList& fc = cyclicGgiPatch_.faceCells();

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


// Return matrix product for coupled boundary
template<class Type>
void cyclicGgiFvPatchField<Type>::updateInterfaceMatrix
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
void cyclicGgiFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const Field<Type>& psiInternal,
    Field<Type>& result,
    const BlockLduMatrix<Type>& m,
    const CoeffField<Type>& coeffs,
    const Pstream::commsTypes commsType
) const
{
    notImplemented
    (
        "void cyclicGgiFvPatchField::initInterfaceMatrixUpdate"
        " for block-coupled matrices of tensor types"
    )
}


template<class Type>
void cyclicGgiFvPatchField<Type>::updateInterfaceMatrix
(
    const Field<Type>&,
    Field<Type>&,
    const BlockLduMatrix<Type>&,
    const CoeffField<Type>&,
    const Pstream::commsTypes commsType
) const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
