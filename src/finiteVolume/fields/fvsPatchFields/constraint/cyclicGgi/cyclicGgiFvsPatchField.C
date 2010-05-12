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

Author
    Martin Beaudoin, Hydro-Quebec, (2008)

Contributor:
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "cyclicGgiFvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
cyclicGgiFvsPatchField<Type>::cyclicGgiFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(p, iF),
    cyclicGgiPatch_(refCast<const cyclicGgiFvPatch>(p))
{}


template<class Type>
cyclicGgiFvsPatchField<Type>::cyclicGgiFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    coupledFvsPatchField<Type>(p, iF, dict, true),
    cyclicGgiPatch_(refCast<const cyclicGgiFvPatch>(p))
{
    if (!isType<cyclicGgiFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "cyclicGgiFvsPatchField<Type>::cyclicGgiFvsPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, surfaceMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not cyclicGgi type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
cyclicGgiFvsPatchField<Type>::cyclicGgiFvsPatchField
(
    const cyclicGgiFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvsPatchField<Type>(ptf, p, iF, mapper),
    cyclicGgiPatch_(refCast<const cyclicGgiFvPatch>(p))
{
    if (!isType<cyclicGgiFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "cyclicGgiFvsPatchField<Type>::cyclicGgiFvsPatchField\n"
            "(\n"
            "    const cyclicGgiFvsPatchField<Type>& ptf,\n"
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
cyclicGgiFvsPatchField<Type>::cyclicGgiFvsPatchField
(
    const cyclicGgiFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(ptf, iF),
    cyclicGgiPatch_(refCast<const cyclicGgiFvPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
