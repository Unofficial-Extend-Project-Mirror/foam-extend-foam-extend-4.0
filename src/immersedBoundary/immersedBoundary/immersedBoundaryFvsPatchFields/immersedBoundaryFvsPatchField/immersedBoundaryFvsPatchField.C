/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "immersedBoundaryFvsPatchField.H"
#include "fvPatchFieldMapper.H"
#include "surfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool immersedBoundaryFvsPatchField<Type>::updateSize()
{
    if (this->patch().size() != this->size())
    {
        this->setSize(this->patch().size());

        // Size has changed
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
immersedBoundaryFvsPatchField<Type>::immersedBoundaryFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    fvsPatchField<Type>(p, iF),
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p))
{}


template<class Type>
immersedBoundaryFvsPatchField<Type>::immersedBoundaryFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    fvsPatchField<Type>(p, iF),   // Do not read base data
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p))
{
    operator=(pTraits<Type>::zero);
}


template<class Type>
immersedBoundaryFvsPatchField<Type>::immersedBoundaryFvsPatchField
(
    const immersedBoundaryFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvsPatchField<Type>(p, iF),  // Do not map base data
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p))
{}


template<class Type>
immersedBoundaryFvsPatchField<Type>::immersedBoundaryFvsPatchField
(
    const immersedBoundaryFvsPatchField<Type>& ptf
)
:
    fvsPatchField<Type>(ptf),
    ibPatch_(ptf.ibPatch())
{}


template<class Type>
immersedBoundaryFvsPatchField<Type>::immersedBoundaryFvsPatchField
(
    const immersedBoundaryFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    fvsPatchField<Type>(ptf, iF),
    ibPatch_(ptf.ibPatch())
{}


template<class Type>
void immersedBoundaryFvsPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    Field<Type>::operator=
    (
        Field<Type>(this->patch().size(), pTraits<Type>::zero)
    );
}


template<class Type>
void immersedBoundaryFvsPatchField<Type>::rmap
(
    const fvsPatchField<Type>& ptf,
    const labelList& addr
)
{
    Field<Type>::operator=
    (
        Field<Type>(this->patch().size(), pTraits<Type>::zero)
    );
}


template<class Type>
void immersedBoundaryFvsPatchField<Type>::updateOnMotion()
{
    if (this->updateSize())
    {
        // Size changed: must reset values.  HJ, 28/Mar/2019
        Field<Type>::operator=
        (
            Field<Type>(this->patch().size(), pTraits<Type>::zero)
        );
    }
}


template<class Type>
void immersedBoundaryFvsPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if (this->updateSize())
    {
        // Size changed: must reset values.  HJ, 28/Mar/2019
        Field<Type>::operator=
        (
            Field<Type>(this->patch().size(), pTraits<Type>::zero)
        );
    }
}


template<class Type>
void immersedBoundaryFvsPatchField<Type>::write(Ostream& os) const
{
    fvsPatchField<Type>::write(os);
    // The value entry needs to be written with zero size
    Field<Type>::null().writeEntry("value", os);
    // this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void immersedBoundaryFvsPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    this->updateSize();
    Field<Type>::operator=(ul);
}


template<class Type>
void immersedBoundaryFvsPatchField<Type>::operator=
(
    const fvsPatchField<Type>& ptf
)
{
    this->check(ptf);
    this->updateSize();
    Field<Type>::operator=(ptf);
}


template<class Type>
void immersedBoundaryFvsPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    this->check(ptf);
    this->updateSize();
    fvsPatchField<Type>::operator=(ptf);
}


template<class Type>
void immersedBoundaryFvsPatchField<Type>::operator=
(
    const Type& t
)
{
    this->updateSize();
    Field<Type>::operator=(t);
}


// Force an assignment, overriding fixedValue status
template<class Type>
void immersedBoundaryFvsPatchField<Type>::operator==
(
    const fvsPatchField<Type>& ptf
)
{
    this->updateSize();
    Field<Type>::operator=(ptf);
}


template<class Type>
void immersedBoundaryFvsPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    this->updateSize();
    Field<Type>::operator=(tf);
}


template<class Type>
void immersedBoundaryFvsPatchField<Type>::operator==
(
    const Type& t
)
{
    this->updateSize();
    Field<Type>::operator=(t);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
