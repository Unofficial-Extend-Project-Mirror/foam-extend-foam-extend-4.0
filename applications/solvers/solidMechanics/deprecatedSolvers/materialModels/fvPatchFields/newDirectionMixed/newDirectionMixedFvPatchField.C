/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    newDirectionMixedFvPatchField

Description
    Doubly mixed fixed value-fixed gradient boundary condition
    separated into a normal and a tangential component given a
    direction vector.  The mixture is controlled by two separate
    valueFraction coefficients in the normal and tangential direction.

\*---------------------------------------------------------------------------*/

#include "newDirectionMixedFvPatchField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void newDirectionMixedFvPatchField<Type>::checkNHat()
{
    scalarField magNHat(Foam::mag(nHat_));

    if (min(magNHat) < SMALL)
    {
        FatalErrorIn("void newDirectionMixedFvPatchField<Type>::checkNHat()")
            << "Incorrectly defined normal direction.  mag = "
            << min(magNHat)
            << abort(FatalError);
    }

    magNHat /= mag(magNHat);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
newDirectionMixedFvPatchField<Type>::newDirectionMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    refValue_(p.size()),
    refGrad_(p.size()),
    nHat_(p.size()),
    normalValueFraction_(p.size()),
    tangentialValueFraction_(p.size())
{}


template<class Type>
newDirectionMixedFvPatchField<Type>::newDirectionMixedFvPatchField
(
    const newDirectionMixedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    refValue_(ptf.refValue_, mapper),
    refGrad_(ptf.refGrad_, mapper),
    nHat_(ptf.nHat_, mapper),
    normalValueFraction_(ptf.normalValueFraction_, mapper),
    tangentialValueFraction_(ptf.tangentialValueFraction_, mapper)
{
    this->checkNHat();
}


template<class Type>
newDirectionMixedFvPatchField<Type>::newDirectionMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict),
    refValue_("refValue", dict, p.size()),
    refGrad_("refGradient", dict, p.size()),
    nHat_("nHat", dict, p.size()),
    normalValueFraction_("normalValueFraction", dict, p.size()),
    tangentialValueFraction_("tangentialValueFraction", dict, p.size())
{
    this->checkNHat();
    evaluate();
}


template<class Type>
newDirectionMixedFvPatchField<Type>::newDirectionMixedFvPatchField
(
    const newDirectionMixedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    nHat_(ptf.nHat_),
    normalValueFraction_(ptf.normalValueFraction_),
    tangentialValueFraction_(ptf.tangentialValueFraction_)
{
    this->checkNHat();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
template<class Type>
void newDirectionMixedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    Field<Type>::autoMap(m);
    refValue_.autoMap(m);
    refGrad_.autoMap(m);
    nHat_.autoMap(m);
    normalValueFraction_.autoMap(m);
    tangentialValueFraction_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
template<class Type>
void newDirectionMixedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const newDirectionMixedFvPatchField<Type>& dmptf =
        refCast<const newDirectionMixedFvPatchField<Type> >(ptf);

    refValue_.rmap(dmptf.refValue_, addr);
    refGrad_.rmap(dmptf.refGrad_, addr);
    nHat_.rmap(dmptf.nHat_, addr);
    normalValueFraction_.rmap(dmptf.normalValueFraction_, addr);
    tangentialValueFraction_.rmap(dmptf.tangentialValueFraction_, addr);
}


// Return gradient at boundary
template<class Type>
tmp<Field<Type> > newDirectionMixedFvPatchField<Type>::snGrad() const
{
    Field<Type> pif = this->patchInternalField();

    const scalarField& deltaCoeffs = this->patch().deltaCoeffs();
    const tensorField nn= nHat_*nHat_;

    Field<Type> normalValue =
        normalValueFraction_*transform(nn, refValue_)
      + (1.0 - normalValueFraction_)*transform(nn, pif + refGrad_/deltaCoeffs);

    Field<Type> tangentialValue =
        tangentialValueFraction_*transform(I - nn, refValue_)
      + (1.0 - tangentialValueFraction_)*
        transform(I - nn, pif + refGrad_/deltaCoeffs);

    return (normalValue + tangentialValue - pif)*deltaCoeffs;
}


// Evaluate the field on the patch
template<class Type>
void newDirectionMixedFvPatchField<Type>::evaluate()
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type> pif = this->patchInternalField();

    const scalarField& deltaCoeffs = this->patch().deltaCoeffs();
    const tensorField nn = nHat_*nHat_;

    Field<Type> normalValue =
        normalValueFraction_*transform(nn, refValue_)
      + (1.0 - normalValueFraction_)*transform(nn, pif + refGrad_/deltaCoeffs);

    Field<Type> tangentialValue =
        tangentialValueFraction_*transform(I - nn, refValue_)
      + (1.0 - tangentialValueFraction_)*
        transform(I - nn, pif + refGrad_/deltaCoeffs);

    Field<Type>::operator=(normalValue + tangentialValue);

    fvPatchField<Type>::evaluate();
}


//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the value of this patchField
template<class Type>
tmp<Field<Type> > newDirectionMixedFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    const Field<Type> one(this->size(), pTraits<Type>::one);
    const tensorField nn= nHat_*nHat_;

    return
        transform(nn, one)*(1.0 - normalValueFraction_)
      + transform(I - nn, one)*(1.0 - tangentialValueFraction_);

}

//- Return the matrix source coefficients corresponding to the
//  evaluation of the value of this patchField
template<class Type>
tmp<Field<Type> > newDirectionMixedFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    const scalarField& deltaCoeffs = this->patch().deltaCoeffs();
    const tensorField nn= nHat_*nHat_;

    return
         normalValueFraction_*transform(nn, refValue_)
       + (1.0 - normalValueFraction_)*transform(nn, refGrad_)/deltaCoeffs
       + tangentialValueFraction_*transform(I - nn, refValue_)
       + (1.0 - tangentialValueFraction_)*
         transform(I - nn, refGrad_)/deltaCoeffs;

    // Alternative; allows fiddling internal/boundary split for value coeffs
//     return
//         *this
//       - scale
//         (
//             valueInternalCoeffs(this->patch().weights()),
//             this->patchInternalField()
//         );
}

//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type> >
newDirectionMixedFvPatchField<Type>::gradientInternalCoeffs() const
{
    const scalarField& deltaCoeffs = this->patch().deltaCoeffs();
    const Field<Type> one(this->size(), pTraits<Type>::one);
    const tensorField nn= nHat_*nHat_;

    return
        -transform(nn, one)*normalValueFraction_*deltaCoeffs
       - transform(I - nn, one)*tangentialValueFraction_*deltaCoeffs;
}

//- Return the matrix source coefficients corresponding to the
//  evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type> >
newDirectionMixedFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    const scalarField& deltaCoeffs = this->patch().deltaCoeffs();
    const tensorField nn= nHat_*nHat_;

    return
         normalValueFraction_*deltaCoeffs*transform(nn, refValue_)
      + (1.0 - normalValueFraction_)*transform(nn, refGrad_)
      +  tangentialValueFraction_*deltaCoeffs*transform(I - nn, refValue_)
      + (1.0 - tangentialValueFraction_)*transform(I - nn, refGrad_);

    // Alternative; allows fiddling internal/boundary split for grad coeffs
//     return
//         snGrad()
//       - scale(gradientInternalCoeffs(), this->patchInternalField());
}


// Write
template<class Type>
void newDirectionMixedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    refValue_.writeEntry("refValue", os);
    refGrad_.writeEntry("refGradient", os);
    nHat_.writeEntry("nHat", os);
    normalValueFraction_.writeEntry("normalValueFraction", os);
    tangentialValueFraction_.writeEntry("tangentialValueFraction", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
