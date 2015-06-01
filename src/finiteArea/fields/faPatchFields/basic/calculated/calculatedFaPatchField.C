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

#include "calculatedFaPatchField.H"
#include "faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
const word& faPatchField<Type>::calculatedType()
{
    return calculatedFaPatchField<Type>::typeName;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
calculatedFaPatchField<Type>::calculatedFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(p, iF)
{}


template<class Type>
calculatedFaPatchField<Type>::calculatedFaPatchField
(
    const calculatedFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
calculatedFaPatchField<Type>::calculatedFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    faPatchField<Type>(p, iF, Field<Type>("value", dict, p.size()))
{}


template<class Type>
calculatedFaPatchField<Type>::calculatedFaPatchField
(
    const calculatedFaPatchField<Type>& ptf
)
:
    faPatchField<Type>(ptf)
{}


template<class Type>
calculatedFaPatchField<Type>::calculatedFaPatchField
(
    const calculatedFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(ptf, iF)
{}


template<class Type>
template<class Type2>
tmp<faPatchField<Type> > faPatchField<Type>::NewCalculatedType
(
    const faPatchField<Type2>& pf
)
{
    typename patchConstructorTable::iterator patchTypeCstrIter =
        patchConstructorTablePtr_->find(pf.patch().type());

    if (patchTypeCstrIter != patchConstructorTablePtr_->end())
    {
        return patchTypeCstrIter()
        (
            pf.patch(),
            DimensionedField<Type, areaMesh>::null()
        );
    }
    else
    {
        return tmp<faPatchField<Type> >
        (
            new calculatedFaPatchField<Type>
            (
                pf.patch(),
                DimensionedField<Type, areaMesh>::null()
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > calculatedFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    FatalErrorIn
    (
        "calculatedFaPatchField<Type>::"
        "valueInternalCoeffs(const tmp<scalarField>&) const"
    )   << "\n    "
           "valueInternalCoeffs cannot be called for a calculatedFaPatchField"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->dimensionedInternalField().name()
        << " in file " << this->dimensionedInternalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "default boundary condition."
        << exit(FatalError);

    return *this;
}


template<class Type>
tmp<Field<Type> > calculatedFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    FatalErrorIn
    (
        "calculatedFaPatchField<Type>::"
        "valueBoundaryCoeffs(const tmp<scalarField>&) const"
    )   << "\n    "
           "valueBoundaryCoeffs cannot be called for a calculatedFaPatchField"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->dimensionedInternalField().name()
        << " in file " << this->dimensionedInternalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "default boundary condition."
        << exit(FatalError);

    return *this;
}


template<class Type>
tmp<Field<Type> > calculatedFaPatchField<Type>::gradientInternalCoeffs() const
{
    FatalErrorIn
    (
        "calculatedFaPatchField<Type>::"
        "gradientInternalCoeffs() const"
    )   << "\n    "
           "gradientInternalCoeffs cannot be called for a "
           "calculatedFaPatchField"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->dimensionedInternalField().name()
        << " in file " << this->dimensionedInternalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "default boundary condition."
        << exit(FatalError);

    return *this;
}


template<class Type>
tmp<Field<Type> > calculatedFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    FatalErrorIn
    (
        "calculatedFaPatchField<Type>::"
        "gradientBoundaryCoeffs() const"
    )   << "\n    "
           "gradientBoundaryCoeffs cannot be called for a "
           "calculatedFaPatchField"
        << "\n    on patch " << this->patch().name()
        << " of field " << this->dimensionedInternalField().name()
        << " in file " << this->dimensionedInternalField().objectPath()
        << "\n    You are probably trying to solve for a field with a "
           "default boundary condition."
        << exit(FatalError);

    return *this;
}


// Write
template<class Type>
void calculatedFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
