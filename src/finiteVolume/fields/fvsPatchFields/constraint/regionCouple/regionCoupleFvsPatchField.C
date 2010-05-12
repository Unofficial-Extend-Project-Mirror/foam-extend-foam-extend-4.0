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

#include "regionCoupleFvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
regionCoupleFvsPatchField<Type>::regionCoupleFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(p, iF),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(p)),
    remoteFieldName_(),
    matrixUpdateBuffer_()
{}


template<class Type>
regionCoupleFvsPatchField<Type>::regionCoupleFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    coupledFvsPatchField<Type>(p, iF, dict),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(p)),
    remoteFieldName_(dict.lookup("remoteField")),
    matrixUpdateBuffer_()
{
    if (!isType<regionCoupleFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "regionCoupleFvsPatchField<Type>::regionCoupleFvsPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, surfaceMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not regionCouple type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
regionCoupleFvsPatchField<Type>::regionCoupleFvsPatchField
(
    const regionCoupleFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvsPatchField<Type>(ptf, p, iF, mapper),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(p)),
    remoteFieldName_(ptf.remoteFieldName_),
    matrixUpdateBuffer_()
{
    if (!isType<regionCoupleFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "regionCoupleFvsPatchField<Type>::regionCoupleFvsPatchField\n"
            "(\n"
            "    const regionCoupleFvsPatchField<Type>& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, surfaceMesh>& iF,\n"
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
regionCoupleFvsPatchField<Type>::regionCoupleFvsPatchField
(
    const regionCoupleFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(ptf, iF),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(ptf.patch())),
    remoteFieldName_(ptf.remoteFieldName_),
    matrixUpdateBuffer_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return neighbour field
template<class Type>
const regionCoupleFvsPatchField<Type>&
regionCoupleFvsPatchField<Type>::shadowPatchField() const
{
    // Lookup neighbour field
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> GeoField;

    const GeoField& coupleField =
        regionCouplePatch_.shadowRegion().
        objectRegistry::lookupObject<GeoField>(remoteFieldName_);

    return refCast<const regionCoupleFvsPatchField<Type> >
    (
        coupleField.boundaryField()[regionCouplePatch_.shadowIndex()]
    );
}


// Write
template<class Type>
void regionCoupleFvsPatchField<Type>::write(Ostream& os) const
{
    fvsPatchField<Type>::write(os);
    os.writeKeyword("remoteField")
        << remoteFieldName_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
