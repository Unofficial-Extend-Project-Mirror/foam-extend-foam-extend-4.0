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

\*---------------------------------------------------------------------------*/

#include "regionCoupleFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
regionCoupleFvPatchField<Type>::regionCoupleFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(p)),
    remoteFieldName_(iF.name()),
    matrixUpdateBuffer_()
{}


template<class Type>
regionCoupleFvPatchField<Type>::regionCoupleFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(p)),
    remoteFieldName_(dict.lookup("remoteField")),
    matrixUpdateBuffer_()
{
    if (!isType<regionCoupleFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "regionCoupleFvPatchField<Type>::regionCoupleFvPatchField\n"
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
}


template<class Type>
regionCoupleFvPatchField<Type>::regionCoupleFvPatchField
(
    const regionCoupleFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(p)),
    remoteFieldName_(ptf.remoteFieldName_),
    matrixUpdateBuffer_()
{
    if (!isType<regionCoupleFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "regionCoupleFvPatchField<Type>::regionCoupleFvPatchField\n"
            "(\n"
            "    const regionCoupleFvPatchField<Type>& ptf,\n"
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
regionCoupleFvPatchField<Type>::regionCoupleFvPatchField
(
    const regionCoupleFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    regionCoupleLduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(ptf.patch())),
    remoteFieldName_(ptf.remoteFieldName_),
    matrixUpdateBuffer_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return neighbour field
template<class Type>
const regionCoupleFvPatchField<Type>&
regionCoupleFvPatchField<Type>::shadowPatchField() const
{
    // Lookup neighbour field
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    const GeoField& coupleField =
        regionCouplePatch_.shadowRegion().
        objectRegistry::lookupObject<GeoField>(remoteFieldName_);

    return refCast<const regionCoupleFvPatchField<Type> >
    (
        coupleField.boundaryField()[regionCouplePatch_.shadowIndex()]
    );
}


// Return neighbour field
template<class Type>
tmp<Field<Type> > regionCoupleFvPatchField<Type>::patchNeighbourField() const
{
    return regionCouplePatch_.interpolate
    (
        shadowPatchField().patchInternalField()
    );
}


template<class Type>
void regionCoupleFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    // Implement weights-based stabilised harmonic interpolation using
    // magnitude of type
    // Algorithm:
    // 1) calculate magnitude of internal field and neighbour field
    // 2) calculate harmonic mean magnitude
    // 3) express harmonic mean magnitude as: mean = w*mOwn + (1 - w)*mNei
    // 4) Based on above, calculate w = (mean - mNei)/(mOwn - mNei)
    // 5) Use weights to interpolate values

    Field<Type> fOwn = this->patchInternalField();
    Field<Type> fNei = this->patchNeighbourField();

    scalarField magFOwn = mag(fOwn);
    scalarField magFNei = mag(fNei);

    // Calculate internal weights using field magnitude
    scalarField weights(fOwn.size());

    forAll (weights, faceI)
    {
        scalar mOwn = magFOwn[faceI];
        scalar mNei = magFNei[faceI];

        scalar den = mOwn - mNei;

        if (mag(den) > SMALL)
        {
            scalar mean = 2.0*mOwn*mNei/(mOwn + mNei);
            weights[faceI] = (mean - mNei)/den;
        }
        else
        {
            weights[faceI] = 0.5;
        }
    }

    // Do interpolation
    Field<Type>::operator=(weights*fOwn + (1.0 - weights)*fNei);
}


// Initialise neighbour processor internal cell data
template<class Type>
void regionCoupleFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField&,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const Pstream::commsTypes
) const
{
    matrixUpdateBuffer_ = this->patch().patchInternalField(psiInternal);
}


// Return matrix product for coupled boundary
template<class Type>
void regionCoupleFvPatchField<Type>::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    scalarField pnf =
        regionCouplePatch_.interpolate
        (
            this->shadowPatchField().matrixUpdateBuffer()
        );

    // Multiply the field by coefficients and add into the result
    const unallocLabelList& fc = regionCouplePatch_.faceCells();

    forAll(fc, elemI)
    {
        result[fc[elemI]] -= coeffs[elemI]*pnf[elemI];
    }
}


// Write
template<class Type>
void regionCoupleFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("remoteField")
        << remoteFieldName_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
