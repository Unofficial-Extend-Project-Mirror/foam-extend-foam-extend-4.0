/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "basicSymmetryFvPatchField.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
basicSymmetryFvPatchField<Type>::basicSymmetryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(p, iF)
{}


template<class Type>
basicSymmetryFvPatchField<Type>::basicSymmetryFvPatchField
(
    const basicSymmetryFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    transformFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
basicSymmetryFvPatchField<Type>::basicSymmetryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    transformFvPatchField<Type>(p, iF, dict)
{
    this->evaluate();
}


template<class Type>
basicSymmetryFvPatchField<Type>::basicSymmetryFvPatchField
(
    const basicSymmetryFvPatchField<Type>& ptf
)
:
    transformFvPatchField<Type>(ptf)
{}


template<class Type>
basicSymmetryFvPatchField<Type>::basicSymmetryFvPatchField
(
    const basicSymmetryFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// return gradient at boundary
template<class Type>
tmp<Field<Type> > basicSymmetryFvPatchField<Type>::snGrad() const
{
    vectorField nHat = this->patch().nf();
    return
    (
        transform(I - 2.0*sqr(nHat), this->patchInternalField())
      - this->patchInternalField()
    )*(this->patch().deltaCoeffs()/2.0);
}


// Evaluate the field on the patch
template<class Type>
void basicSymmetryFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    vectorField nHat = this->patch().nf();

    Field<Type> pif = this->patchInternalField();
    Field<Type>::operator=
    (
        0.5*(pif + transform(I - 2.0*sqr(nHat), pif))
    );

    transformFvPatchField<Type>::evaluate();
}


// Return defining fields
template<class Type>
tmp<Field<Type> > basicSymmetryFvPatchField<Type>::snGradTransformDiag() const
{
    vectorField nHat = this->patch().nf();
    vectorField diag(nHat.size());

    diag.replace(vector::X, mag(nHat.component(vector::X)));
    diag.replace(vector::Y, mag(nHat.component(vector::Y)));
    diag.replace(vector::Z, mag(nHat.component(vector::Z)));

    return transformFieldMask<Type>(pow<vector, pTraits<Type>::rank>(diag));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
