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

\*---------------------------------------------------------------------------*/

#include "regionCouplingFvPatchField.H"
#include "symmTransformField.H"
#include "harmonic.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
const Foam::Field<Type>&
regionCouplingFvPatchField<Type>::originalPatchField() const
{
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        // Store original field for symmetric evaluation
        // Henrik Rusche, Aug/2011

        originalPatchField_ = *this;
        curTimeIndex_ = this->db().time().timeIndex();
    }

    return originalPatchField_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
regionCouplingFvPatchField<Type>::regionCouplingFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(p)),
    remoteFieldName_(iF.name()),
    matrixUpdateBuffer_(),
    originalPatchField_(),
    curTimeIndex_(-1)
{}


template<class Type>
regionCouplingFvPatchField<Type>::regionCouplingFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(p)),
    remoteFieldName_(dict.lookupOrDefault<word>("remoteField", iF.name())),
    matrixUpdateBuffer_(),
    originalPatchField_(),
    curTimeIndex_(-1)
{
    if (!isType<regionCoupleFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "regionCouplingFvPatchField<Type>::regionCouplingFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not regionCouple type. "
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
regionCouplingFvPatchField<Type>::regionCouplingFvPatchField
(
    const regionCouplingFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(p)),
    remoteFieldName_(ptf.remoteFieldName_),
    matrixUpdateBuffer_(),
    originalPatchField_(),
    curTimeIndex_(-1)
{
    if (!isType<regionCoupleFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "regionCouplingFvPatchField<Type>::regionCouplingFvPatchField\n"
            "(\n"
            "    const regionCouplingFvPatchField<Type>& ptf,\n"
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
regionCouplingFvPatchField<Type>::regionCouplingFvPatchField
(
    const regionCouplingFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    ggiLduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(ptf.patch())),
    remoteFieldName_(ptf.remoteFieldName_),
    matrixUpdateBuffer_(),
    originalPatchField_(),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return a named shadow patch field
template<class Type>
template<class LookupField, class LookupType>
const typename LookupField::PatchFieldType&
regionCouplingFvPatchField<Type>::lookupShadowPatchField
(
    const word& name,
    const LookupField*,
    const LookupType*
) const
{
    // Lookup neighbour field
    const LookupField& shadowField =
        regionCouplePatch_.shadowRegion().
        objectRegistry::template lookupObject<LookupField>(name);

    return shadowField.boundaryField()[regionCouplePatch_.shadowIndex()];
}


// Return shadow patch field
template<class Type>
const regionCouplingFvPatchField<Type>&
regionCouplingFvPatchField<Type>::shadowPatchField() const
{
    // Lookup neighbour field
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    return refCast<const regionCouplingFvPatchField<Type> >
    (
        lookupShadowPatchField<GeoField, Type>(remoteFieldName_)
    );
}


// Return neighbour field
template<class Type>
tmp<Field<Type> > regionCouplingFvPatchField<Type>::patchNeighbourField() const
{
    Field<Type> sField = shadowPatchField().patchInternalField();

    tmp<Field<Type> > tpnf
    (
         regionCouplePatch_.interpolate
         (
             shadowPatchField().patchInternalField()
         )
    );

    Field<Type>& pnf = tpnf();

    if (regionCouplePatch_.bridgeOverlap())
    {
        // Symmetry treatment used for overlap
        vectorField nHat = this->patch().nf();

        // Use mirrored neighbour field for interpolation
        // HJ, 21/Jan/2009
        Field<Type> bridgeField =
            transform(I - 2.0*sqr(nHat), this->patchInternalField());

        regionCouplePatch_.bridge(bridgeField, pnf);
    }

    return tpnf;
}


// Return neighbour field given internal cell data
template<class Type>
tmp<Field<Type> > regionCouplingFvPatchField<Type>::patchNeighbourField
(
    const word& name
) const
{
    // Lookup neighbour field
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    return regionCouplePatch_.interpolate
    (
        lookupShadowPatchField<GeoField, Type>(name).patchInternalField()
    );

    // Note: this field is not bridged because local data does not exist
    // for named field.  HJ, 27/Sep/2011
}


template<class Type>
void regionCouplingFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if(debug)
    {
        Info << "In regionCouplingFvPatchField<Type>::initEvaluate() on "
            << this->dimensionedInternalField().name()
            << " in " << this->patch().boundaryMesh().mesh().name()
            << " " << this->updated() << endl;
    }

    // Interpolation must happen at init

    // Note: If used with interpolation - either on explicitly or called by the
    // laplacian operator, the values set here are overridden by the
    // interpolation scheme. In order to get the same diffusivities on
    //  both sides an identical interpolation scheme must be used.
    // Note^2: Even if harmonic used, the interpolation is still wrong for most
    // CHT cases since (cell values vs. face values)
    // Note^3: None of this is intuitiv - fix requires low-level changes!
    // HR, 8/Jun/2012

    const Field<Type>& fOwn = this->originalPatchField();
    const Field<Type> fNei = regionCouplePatch_.interpolate
    (
        this->shadowPatchField().originalPatchField()
    );

    // Do interpolation
    harmonic<Type> interp(this->patch().boundaryMesh().mesh());

    scalarField weights = interp.weights(fOwn, fNei, this->patch());

    Field<Type>::operator=(weights*fOwn + (1.0 - weights)*fNei);

    if (regionCouplePatch_.bridgeOverlap())
    {
        // Symmetry treatment used for overlap
        vectorField nHat = this->patch().nf();

        Field<Type> pif = this->patchInternalField();

        Field<Type> bridgeField =
            0.5*(pif + transform(I - 2.0*sqr(nHat), pif));

        regionCouplePatch_.bridge(bridgeField, *this);
    }
}


template<class Type>
void regionCouplingFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    // No interpolation allowed

    fvPatchField<Type>::evaluate();
}


template<class Type>
void regionCouplingFvPatchField<Type>::updateCoeffs()
{
    if(debug)
    {
        Info << "In regionCouplingFvPatchField<Type>::updateCoeffs() on "
            << this->dimensionedInternalField().name()
            << " in " << this->patch().boundaryMesh().mesh().name()
            << " " << this->updated() << endl;
    }

    if (this->updated())
    {
        return;
    }

    Field<Type> fOwn = this->patchInternalField();
    Field<Type> fNei = this->patchNeighbourField();

    // Do interpolation
    harmonic<Type> interp(this->patch().boundaryMesh().mesh());

    scalarField weights = interp.weights(fOwn, fNei, this->patch());

    Field<Type>::operator=(weights*fOwn + (1.0 - weights)*fNei);

    if (regionCouplePatch_.bridgeOverlap())
    {
        // Symmetry treatment used for overlap
        vectorField nHat = this->patch().nf();

        Field<Type> pif = this->patchInternalField();

        Field<Type> bridgeField =
            0.5*(pif + transform(I - 2.0*sqr(nHat), pif));

        regionCouplePatch_.bridge(bridgeField, *this);
    }
}


template<class Type>
tmp<Field<Type> > regionCouplingFvPatchField<Type>::snGrad() const
{
    if (regionCouplePatch_.coupled())
    {
        return coupledFvPatchField<Type>::snGrad();
    }
    else
    {
        return fvPatchField<Type>::snGrad();
    }
}


// Initialise neighbour processor internal cell data
template<class Type>
void regionCouplingFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes,
    const bool switchToLhs
) const
{
    if (regionCouplePatch_.coupled())
    {
        // Prepare local matrix update buffer for the remote side.
        // Note that only remote side will have access to its psiInternal
        // as they are on different regions

        // Since interpolation needs to happen on the shadow, and within the
        // init, prepare interpolation for the other side.
        matrixUpdateBuffer_ =
            this->shadowPatchField().regionCouplePatch().interpolate
            (
                this->patch().patchInternalField(psiInternal)
            );
    }
    else
    {
        FatalErrorIn
        (
            "regionCouplingFvPatchField<Type>::initInterfaceMatrixUpdate"
        )   << "init matrix update called in detached state"
            << abort(FatalError);
    }
}


// Return matrix product for coupled boundary
template<class Type>
void regionCouplingFvPatchField<Type>::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction ,
    const Pstream::commsTypes,
    const bool switchToLhs
) const
{
    if (regionCouplePatch_.coupled())
    {
        // Note: interpolation involves parallel communications and needs to
        // happen during init.  This changes the use of matrix update buffer
        // compared to earlier versions
        // HJ, 28/Sep/2011
        scalarField pnf = this->shadowPatchField().matrixUpdateBuffer();

        // Multiply the field by coefficients and add into the result
        const unallocLabelList& fc = regionCouplePatch_.faceCells();

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
    else
    {
        FatalErrorIn
        (
            "regionCouplingFvPatchField<Type>::updateInterfaceMatrix"
        )   << "Matrix update called in detached state"
            << abort(FatalError);
    }
}


// Write
template<class Type>
void regionCouplingFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("remoteField")
        << remoteFieldName_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
