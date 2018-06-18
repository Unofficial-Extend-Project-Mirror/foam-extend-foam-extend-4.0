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

#include "oversetFvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
oversetFvsPatchField<Type>::oversetFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(p, iF),
    oversetPatch_(refCast<const oversetFvPatch>(p))
{}


template<class Type>
oversetFvsPatchField<Type>::oversetFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const Field<Type>& f
)
:
    coupledFvsPatchField<Type>(p, iF, f),
    oversetPatch_(refCast<const oversetFvPatch>(p))
{}


// Construct by mapping given oversetFvsPatchField<Type>
template<class Type>
oversetFvsPatchField<Type>::oversetFvsPatchField
(
    const oversetFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvsPatchField<Type>(ptf, p, iF, mapper),
    oversetPatch_(refCast<const oversetFvPatch>(p))
{
    if (!isType<oversetFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "oversetFvsPatchField<Type>::oversetFvsPatchField\n"
            "(\n"
            "    const oversetFvsPatchField<Type>& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, surfaceMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
oversetFvsPatchField<Type>::oversetFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    coupledFvsPatchField<Type>(p, iF, dict),
    oversetPatch_(refCast<const oversetFvPatch>(p))
{
    if (!isType<oversetFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "oversetFvsPatchField<Type>::oversetFvsPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not overset type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
oversetFvsPatchField<Type>::oversetFvsPatchField
(
    const oversetFvsPatchField<Type>& ptf
)
:
    coupledFvsPatchField<Type>(ptf),
    oversetPatch_(refCast<const oversetFvPatch>(ptf.patch()))
{}


template<class Type>
oversetFvsPatchField<Type>::oversetFvsPatchField
(
    const oversetFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(ptf, iF),
    oversetPatch_(refCast<const oversetFvPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
oversetFvsPatchField<Type>::~oversetFvsPatchField()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
