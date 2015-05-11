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

#include "immersedBoundaryFvsPatchField.H"
#include "fvPatchFieldMapper.H"
#include "surfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


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
    fvsPatchField<Type>(p, iF, dict),
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p))
{}


template<class Type>
immersedBoundaryFvsPatchField<Type>::immersedBoundaryFvsPatchField
(
    const immersedBoundaryFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvsPatchField<Type>(ptf, p, iF, mapper),
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


// template<class Type>
// void immersedBoundaryFvsPatchField<Type>::operator=
// (
//     const fvPatchField<Type>& ptf
// )
// {
//     const immersedBoundaryFvPatchField<Type>& ibf =
//         refCast<const immersedBoundaryFvPatchField<Type> > (ptf);

//     this->check(ptf);
//     fvsPatchField<Type>::operator=(ptf);
// }


template<class Type>
void immersedBoundaryFvsPatchField<Type>::write(Ostream& os) const
{
    fvsPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
