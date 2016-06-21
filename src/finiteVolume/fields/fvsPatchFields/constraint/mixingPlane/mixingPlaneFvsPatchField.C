/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

Description
    Mixing plane fvs patch

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "mixingPlaneFvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
mixingPlaneFvsPatchField<Type>::mixingPlaneFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(p, iF),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(p))
{}


template<class Type>
mixingPlaneFvsPatchField<Type>::mixingPlaneFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    coupledFvsPatchField<Type>(p, iF, dict, true),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(p))
{
    if (!isType<mixingPlaneFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "mixingPlaneFvsPatchField<Type>::mixingPlaneFvsPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, surfaceMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not mixingPlane type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
mixingPlaneFvsPatchField<Type>::mixingPlaneFvsPatchField
(
    const mixingPlaneFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvsPatchField<Type>(ptf, p, iF, mapper),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(p))
{
    if (!isType<mixingPlaneFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "mixingPlaneFvsPatchField<Type>::mixingPlaneFvsPatchField\n"
            "(\n"
            "    const mixingPlaneFvsPatchField<Type>& ptf,\n"
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
mixingPlaneFvsPatchField<Type>::mixingPlaneFvsPatchField
(
    const mixingPlaneFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(ptf, iF),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
