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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description

\*---------------------------------------------------------------------------*/

#include "fixedGradientCorrectedFvPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
fixedGradientCorrectedFvPatchField<Type>::fixedGradientCorrectedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    correctedFvPatchField<Type>(p, iF),
    gradient_(p.size(), pTraits<Type>::zero)
{}


template<class Type>
fixedGradientCorrectedFvPatchField<Type>::fixedGradientCorrectedFvPatchField
(
    const fixedGradientCorrectedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    correctedFvPatchField<Type>(ptf, p, iF, mapper),
    gradient_(ptf.gradient_, mapper)
{}


template<class Type>
fixedGradientCorrectedFvPatchField<Type>::fixedGradientCorrectedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    correctedFvPatchField<Type>(p, iF, dict),
    gradient_("gradient", dict, p.size())
{
    if (dict.found("value"))
    {
        Field<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        Field<Type>::operator=
        (
            this->patchInternalField() 
          + gradient_/this->patch().deltaCoeffs()
          + this->corrVecGrad()
        );
    }

    this->nGradInternal() = gradient_;
}


template<class Type>
fixedGradientCorrectedFvPatchField<Type>::fixedGradientCorrectedFvPatchField
(
    const fixedGradientCorrectedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    correctedFvPatchField<Type>(ptf, iF),
    gradient_(ptf.gradient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
template<class Type>
void fixedGradientCorrectedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    Field<Type>::autoMap(m);
    gradient_.autoMap(m);

//     if (m.resizeOnly())
//     {
//         setSize(m.size());
//         gradient_.setSize(m.size());
//     }
//     else
//     {
//         Field<Type>::autoMap((const FieldMapper&)m);
//         gradient_.autoMap((const FieldMapper&)m);
//     }
}


// Reverse-map the given fvPatchField onto this fvPatchField
template<class Type>
void fixedGradientCorrectedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const fixedGradientCorrectedFvPatchField<Type>& fgptf =
        refCast<const fixedGradientCorrectedFvPatchField<Type> >(ptf);

    gradient_.rmap(fgptf.gradient_, addr);
}


// Evaluate the field on the patch
template<class Type>
void fixedGradientCorrectedFvPatchField<Type>::evaluate()
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        this->patchInternalField()
      + this->corrVecGrad()
      + gradient_/this->patch().deltaCoeffs()
    );

//     Field<Type>::operator=
//     (
//         this->patchInternalField()
//       + this->corrVecGrad()
//       + 0.5*gradient_/this->patch().deltaCoeffs()
//       + 0.5*this->nGradInternal()/this->patch().deltaCoeffs()
//     );

//     Field<Type> newValue
//     (
//         this->patchInternalField()
//       + this->corrVecGrad()
//       + 0.5*gradient_/this->patch().deltaCoeffs()
//       + 0.5*this->nGradInternal()/this->patch().deltaCoeffs()
//     );

//     scalar alpha = 0.25;

//     newValue = alpha*newValue + (1.0 - alpha)*(*this);

//     Field<Type>::operator=(newValue);

    fvPatchField<Type>::evaluate();
}


//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the value of this patchField
template<class Type>
tmp<Field<Type> > fixedGradientCorrectedFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type> >(new Field<Type>(this->size(), pTraits<Type>::one));
}

//- Return the matrix source coefficients corresponding to the
//  evaluation of the value of this patchField
template<class Type>
tmp<Field<Type> > fixedGradientCorrectedFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
//     return this->corrVecGrad()
//       + 0.5*gradient()/this->patch().deltaCoeffs()
//       + 0.5*this->nGradInternal()/this->patch().deltaCoeffs();

    return this->corrVecGrad()
        + gradient()/this->patch().deltaCoeffs();
}

//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type> > fixedGradientCorrectedFvPatchField<Type>::
gradientInternalCoeffs() const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}

//- Return the matrix source coefficients corresponding to the
//  evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type> > fixedGradientCorrectedFvPatchField<Type>::
gradientBoundaryCoeffs() const
{
    return gradient();
}


// Write
template<class Type>
void fixedGradientCorrectedFvPatchField<Type>::write(Ostream& os) const
{
    correctedFvPatchField<Type>::write(os);
    gradient_.writeEntry("gradient", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
