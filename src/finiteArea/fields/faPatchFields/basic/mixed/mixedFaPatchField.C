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

\*---------------------------------------------------------------------------*/

#include "mixedFaPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
mixedFaPatchField<Type>::mixedFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(p, iF),
    refValue_(p.size()),
    refGrad_(p.size()),
    valueFraction_(p.size())
{}


template<class Type>
mixedFaPatchField<Type>::mixedFaPatchField
(
    const mixedFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faPatchField<Type>(ptf, p, iF, mapper),
    refValue_(ptf.refValue_, mapper),
    refGrad_(ptf.refGrad_, mapper),
    valueFraction_(ptf.valueFraction_, mapper)
{}


template<class Type>
mixedFaPatchField<Type>::mixedFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    faPatchField<Type>(p, iF),
    refValue_("refValue", dict, p.size()),
    refGrad_("refGradient", dict, p.size()),
    valueFraction_("valueFraction", dict, p.size())
{
    evaluate();
}


template<class Type>
mixedFaPatchField<Type>::mixedFaPatchField
(
    const mixedFaPatchField<Type>& ptf
)
:
    faPatchField<Type>(ptf),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_)
{}


template<class Type>
mixedFaPatchField<Type>::mixedFaPatchField
(
    const mixedFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(ptf, iF),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void mixedFaPatchField<Type>::autoMap
(
    const faPatchFieldMapper& m
)
{
    Field<Type>::autoMap(m);
    refValue_.autoMap(m);
    refGrad_.autoMap(m);
    valueFraction_.autoMap(m);
}


template<class Type>
void mixedFaPatchField<Type>::rmap
(
    const faPatchField<Type>& ptf,
    const labelList& addr
)
{
    faPatchField<Type>::rmap(ptf, addr);

    const mixedFaPatchField<Type>& mptf =
        refCast<const mixedFaPatchField<Type> >(ptf);

    refValue_.rmap(mptf.refValue_, addr);
    refGrad_.rmap(mptf.refGrad_, addr);
    valueFraction_.rmap(mptf.valueFraction_, addr);
}


template<class Type>
void mixedFaPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        valueFraction_*refValue_
      +
        (1.0 - valueFraction_)*
        (
            this->patchInternalField()
          + refGrad_/this->patch().deltaCoeffs()
        )
    );

    faPatchField<Type>::evaluate();
}


template<class Type>
tmp<Field<Type> > mixedFaPatchField<Type>::snGrad() const
{
    return
        valueFraction_
       *(refValue_ - this->patchInternalField())
       *this->patch().deltaCoeffs()
      + (1.0 - valueFraction_)*refGrad_;
}


template<class Type>
tmp<Field<Type> > mixedFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return Type(pTraits<Type>::one)*(1.0 - valueFraction_);

}


template<class Type>
tmp<Field<Type> > mixedFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return
         valueFraction_*refValue_
       + (1.0 - valueFraction_)*refGrad_/this->patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type> > mixedFaPatchField<Type>::gradientInternalCoeffs() const
{
    return -Type(pTraits<Type>::one)*valueFraction_*this->patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type> > mixedFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    return
        valueFraction_*this->patch().deltaCoeffs()*refValue_
      + (1.0 - valueFraction_)*refGrad_;
}


template<class Type>
void mixedFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    refValue_.writeEntry("refValue", os);
    refGrad_.writeEntry("refGradient", os);
    valueFraction_.writeEntry("valueFraction", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
