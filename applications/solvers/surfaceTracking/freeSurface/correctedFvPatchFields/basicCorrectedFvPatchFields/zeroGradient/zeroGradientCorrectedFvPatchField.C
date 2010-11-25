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

#include "zeroGradientCorrectedFvPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
zeroGradientCorrectedFvPatchField<Type>::zeroGradientCorrectedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    correctedFvPatchField<Type>(p, iF)
{}


template<class Type>
zeroGradientCorrectedFvPatchField<Type>::zeroGradientCorrectedFvPatchField
(
    const zeroGradientCorrectedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    correctedFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
zeroGradientCorrectedFvPatchField<Type>::zeroGradientCorrectedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    correctedFvPatchField<Type>(p, iF, dict)
{
    Field<Type>::operator=
    (
        this->patchInternalField() + this->corrVecGrad()
    );
}


template<class Type>
zeroGradientCorrectedFvPatchField<Type>::zeroGradientCorrectedFvPatchField
(
    const zeroGradientCorrectedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    correctedFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
template<class Type>
void zeroGradientCorrectedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    Field<Type>::autoMap(m);

//     if (m.resizeOnly())
//     {
//         setSize(m.size());
//     }
//     else
//     {
//         Field<Type>::autoMap((const FieldMapper&)m);
//     }
}


// Reverse-map the given fvPatchField onto this fvPatchField
template<class Type>
void zeroGradientCorrectedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);
}


// Evaluate the field on the patch
template<class Type>
void zeroGradientCorrectedFvPatchField<Type>::evaluate()
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        this->patchInternalField()
      + this->corrVecGrad()
    );

    fvPatchField<Type>::evaluate();
}


//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the value of this patchField
template<class Type>
tmp<Field<Type> > zeroGradientCorrectedFvPatchField<Type>::
valueInternalCoeffs(const tmp<scalarField>&) const
{
    return tmp<Field<Type> >(new Field<Type>(this->size(), pTraits<Type>::one));
}

//- Return the matrix source coefficients corresponding to the
//  evaluation of the value of this patchField
template<class Type>
tmp<Field<Type> > zeroGradientCorrectedFvPatchField<Type>::
valueBoundaryCoeffs(const tmp<scalarField>&) const
{
    return tmp<Field<Type> >
    (
//         new Field<Type>(this->size(), pTraits<Type>::zero)
        new Field<Type>(this->corrVecGrad())
    );
}

//- Return the matrix diagonal coefficients corresponding to the
//  evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type> > zeroGradientCorrectedFvPatchField<Type>::
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
tmp<Field<Type> > zeroGradientCorrectedFvPatchField<Type>::
gradientBoundaryCoeffs() const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}


// Write
template<class Type>
void zeroGradientCorrectedFvPatchField<Type>::write(Ostream& os) const
{
    correctedFvPatchField<Type>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
