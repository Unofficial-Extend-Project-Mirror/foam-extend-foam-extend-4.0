/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

\*---------------------------------------------------------------------------*/

#include "coupledFaPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
coupledFaPatchField<Type>::coupledFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    lduInterfaceField(refCast<const lduInterface>(p)),
    faPatchField<Type>(p, iF)
{}


template<class Type>
coupledFaPatchField<Type>::coupledFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const Field<Type>& f
)
:
    lduInterfaceField(refCast<const lduInterface>(p)),
    faPatchField<Type>(p, iF, f)
{}


template<class Type>
coupledFaPatchField<Type>::coupledFaPatchField
(
    const coupledFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    lduInterfaceField(refCast<const lduInterface>(p)),
    faPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
coupledFaPatchField<Type>::coupledFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    lduInterfaceField(refCast<const lduInterface>(p)),
    faPatchField<Type>(p, iF, dict)
{}


template<class Type>
coupledFaPatchField<Type>::coupledFaPatchField
(
    const coupledFaPatchField<Type>& ptf
)
:
    lduInterfaceField(refCast<const lduInterface>(ptf.patch())),
    faPatchField<Type>(ptf)
{}


template<class Type>
coupledFaPatchField<Type>::coupledFaPatchField
(
    const coupledFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    lduInterfaceField(refCast<const lduInterface>(ptf.patch())),
    faPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > coupledFaPatchField<Type>::snGrad() const
{
    return
        (patchNeighbourField() - this->patchInternalField())
       *this->patch().deltaCoeffs();
}


template<class Type>
void coupledFaPatchField<Type>::initEvaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
}


template<class Type>
void coupledFaPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    Field<Type>::operator=
    (
        this->patch().weights()*this->patchInternalField()
      + (1.0 - this->patch().weights())*patchNeighbourField()
    );
}


template<class Type>
tmp<Field<Type> > coupledFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    return Type(pTraits<Type>::one)*w;
}


template<class Type>
tmp<Field<Type> > coupledFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    return Type(pTraits<Type>::one)*(1.0 - w);
}


template<class Type>
tmp<Field<Type> > coupledFaPatchField<Type>::gradientInternalCoeffs() const
{
    return -Type(pTraits<Type>::one)*this->patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type> > coupledFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    return -gradientInternalCoeffs();
}


template<class Type>
void coupledFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
