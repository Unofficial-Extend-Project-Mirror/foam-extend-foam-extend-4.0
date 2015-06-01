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

#include "IOobject.H"
#include "dictionary.H"
#include "faMesh.H"
#include "faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
faePatchField<Type>::faePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF)
{}


template<class Type>
faePatchField<Type>::faePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const Field<Type>& f
)
:
    Field<Type>(f),
    patch_(p),
    internalField_(iF)
{}


template<class Type>
faePatchField<Type>::faePatchField
(
    const faePatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    Field<Type>(ptf, mapper),
    patch_(p),
    internalField_(iF)
{}


template<class Type>
faePatchField<Type>::faePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const dictionary& dict
)
:
    Field<Type>(p.size()),
    patch_(p),
    internalField_(iF)
{
    if (dict.found("value"))
    {
        faePatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        faePatchField<Type>::operator=(pTraits<Type>::zero);
    }
}


template<class Type>
faePatchField<Type>::faePatchField
(
    const faePatchField<Type>& ptf
)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(ptf.internalField_)
{}


template<class Type>
faePatchField<Type>::faePatchField
(
    const faePatchField<Type>& ptf,
    const DimensionedField<Type, edgeMesh>& iF
)
:
    Field<Type>(ptf),
    patch_(ptf.patch_),
    internalField_(iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const objectRegistry& faePatchField<Type>::db() const
{
    //HR 12.3.10: Lookup fields from the field DB rather than the mesh
    return internalField_.db();
}


template<class Type>
template<class GeometricField, class Type2>
const typename GeometricField::PatchFieldType& Foam::faePatchField<Type>::lookupPatchField
(
    const word& name,
    const GeometricField*,
    const Type2*
) const
{
    return patch_.patchField<GeometricField, Type2>
    (
        internalField_.db().objectRegistry::template
        lookupObject<GeometricField>(name)
    );
}


template<class Type>
void faePatchField<Type>::check(const faePatchField<Type>& ptf) const
{
    if (&patch_ != &(ptf.patch_))
    {
        FatalErrorIn("PatchField<Type>::check(const faePatchField<Type>&)")
            << "different patches for faePatchField<Type>s"
            << abort(FatalError);
    }
}


// Map from self
template<class Type>
void faePatchField<Type>::autoMap
(
    const faPatchFieldMapper& m
)
{
    Field<Type>::autoMap(m);
}


// Reverse-map the given faePatchField onto this faePatchField
template<class Type>
void faePatchField<Type>::rmap
(
    const faePatchField<Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::rmap(ptf, addr);
}


// Write
template<class Type>
void faePatchField<Type>::write(Ostream& os) const
{
    os.writeKeyword("type") << type() << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void faePatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    Field<Type>::operator=(ul);
}


template<class Type>
void faePatchField<Type>::operator=
(
    const faePatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator=(ptf);
}


template<class Type>
void faePatchField<Type>::operator+=
(
    const faePatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator+=(ptf);
}


template<class Type>
void faePatchField<Type>::operator-=
(
    const faePatchField<Type>& ptf
)
{
    check(ptf);
    Field<Type>::operator-=(ptf);
}


template<class Type>
void faePatchField<Type>::operator*=
(
    const faePatchField<scalar>& ptf
)
{
    if (&patch_ != &ptf.patch())
    {
        FatalErrorIn
        (
            "PatchField<Type>::operator*=(const faePatchField<scalar>& ptf)"
        )   << "incompatible patches for patch fields"
            << abort(FatalError);
    }

    Field<Type>::operator*=(ptf);
}


template<class Type>
void faePatchField<Type>::operator/=
(
    const faePatchField<scalar>& ptf
)
{
    if (&patch_ != &ptf.patch())
    {
        FatalErrorIn
        (
            "PatchField<Type>::operator/=(const faePatchField<scalar>& ptf)"
        )   << "    incompatible patches for patch fields"
            << abort(FatalError);
    }

    Field<Type>::operator/=(ptf);
}


template<class Type>
void faePatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator+=(tf);
}


template<class Type>
void faePatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator-=(tf);
}


template<class Type>
void faePatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    Field<Type>::operator*=(tf);
}


template<class Type>
void faePatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    Field<Type>::operator/=(tf);
}


template<class Type>
void faePatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


template<class Type>
void faePatchField<Type>::operator+=
(
    const Type& t
)
{
    Field<Type>::operator+=(t);
}


template<class Type>
void faePatchField<Type>::operator-=
(
    const Type& t
)
{
    Field<Type>::operator-=(t);
}


template<class Type>
void faePatchField<Type>::operator*=
(
    const scalar s
)
{
    Field<Type>::operator*=(s);
}


template<class Type>
void faePatchField<Type>::operator/=
(
    const scalar s
)
{
    Field<Type>::operator/=(s);
}


// Force an assignment, overriding fixedValue status
template<class Type>
void faePatchField<Type>::operator==
(
    const faePatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
}


template<class Type>
void faePatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
}


template<class Type>
void faePatchField<Type>::operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Ostream& operator<<(Ostream& os, const faePatchField<Type>& ptf)
{
    ptf.write(os);

    os.check("Ostream& operator<<(Ostream&, const faePatchField<Type>&");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "newFaePatchField.C"

// ************************************************************************* //
