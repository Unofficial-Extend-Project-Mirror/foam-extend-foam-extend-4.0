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

#include "calculatedFaePatchField.H"
#include "faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
const word& faePatchField<Type>::calculatedType()
{
    return calculatedFaePatchField<Type>::typeName;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
calculatedFaePatchField<Type>::calculatedFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF
)
:
    faePatchField<Type>(p, iF)
{}


template<class Type>
calculatedFaePatchField<Type>::calculatedFaePatchField
(
    const calculatedFaePatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faePatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
calculatedFaePatchField<Type>::calculatedFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const dictionary& dict
)
:
    faePatchField<Type>(p, iF, Field<Type>("value", dict, p.size()))
{}


template<class Type>
calculatedFaePatchField<Type>::calculatedFaePatchField
(
    const calculatedFaePatchField<Type>& ptf
)
:
    faePatchField<Type>(ptf)
{}


template<class Type>
calculatedFaePatchField<Type>::calculatedFaePatchField
(
    const calculatedFaePatchField<Type>& ptf,
    const DimensionedField<Type, edgeMesh>& iF
)
:
    faePatchField<Type>(ptf, iF)
{}


template<class Type>
template<class Type2>
tmp<faePatchField<Type> > faePatchField<Type>::NewCalculatedType
(
    const faePatchField<Type2>& pf
)
{
    typename patchConstructorTable::iterator patchTypeCstrIter =
        patchConstructorTablePtr_->find(pf.patch().type());

    if (patchTypeCstrIter != patchConstructorTablePtr_->end())
    {
        return patchTypeCstrIter()
        (
            pf.patch(),
            DimensionedField<Type, edgeMesh>::null()
        );
    }
    else
    {
        return tmp<faePatchField<Type> >
        (
            new calculatedFaePatchField<Type>
            (
                pf.patch(),
                DimensionedField<Type, edgeMesh>::null()
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Write
template<class Type>
void calculatedFaePatchField<Type>::write(Ostream& os) const
{
    faePatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
