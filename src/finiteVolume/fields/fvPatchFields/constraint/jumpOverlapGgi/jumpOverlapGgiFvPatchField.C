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

Author
    Ilaria De Dominicis, General Electric Power, (March 2016)

Contributor
    Hrvoje Jasak, Wikki Ltd.

GE CONFIDENTIAL INFORMATION 2016 General Electric Company. All Rights Reserved

\*---------------------------------------------------------------------------*/

#include "jumpOverlapGgiFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
jumpOverlapGgiFvPatchField<Type>::jumpOverlapGgiFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    overlapGgiFvPatchField<Type>(p, iF)
{}


template<class Type>
jumpOverlapGgiFvPatchField<Type>::jumpOverlapGgiFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    overlapGgiFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
jumpOverlapGgiFvPatchField<Type>::jumpOverlapGgiFvPatchField
(
    const jumpOverlapGgiFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    overlapGgiFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
jumpOverlapGgiFvPatchField<Type>::jumpOverlapGgiFvPatchField
(
    const jumpOverlapGgiFvPatchField<Type>& ptf
)
:
    overlapGgiFvPatchField<Type>(ptf)
{}


template<class Type>
jumpOverlapGgiFvPatchField<Type>::jumpOverlapGgiFvPatchField
(
    const jumpOverlapGgiFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    overlapGgiFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
tmp<Field<Type> >
jumpOverlapGgiFvPatchField<Type>::patchNeighbourField() const
{
    return overlapGgiFvPatchField<Type>::patchNeighbourField() + jump();
}


template<class Type>
void jumpOverlapGgiFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{

    // Get shadow face-cells and assemble shadow field
    const unallocLabelList& sfc = this->overlapGgiPatch().shadow().faceCells();

    scalarField sField(sfc.size());

    forAll (sField, i)
    {
        sField[i] = psiInternal[sfc[i]];
    }

    scalarField pnf = this->overlapGgiPatch().interpolate(sField);

    // Multiply the field by coefficients and add into the result
    const unallocLabelList& fc = this->overlapGgiPatch().faceCells();

    if
    (
        reinterpret_cast<const void*>(&psiInternal)
     == reinterpret_cast<const void*>(&this->internalField())
    )
    {
        const Field<scalar> jf = jump()().component(cmpt);

        if (switchToLhs)
        {
            forAll(fc, elemI)
            {
                 result[fc[elemI]] += coeffs[elemI]*(pnf[elemI] + jf[elemI]);
            }
        }
        else
        {
            forAll(fc, elemI)
            {
                result[fc[elemI]] -= coeffs[elemI]*(pnf[elemI] + jf[elemI]);
            }
        }
    }
    else
    {
        if (switchToLhs)
        {
            forAll(fc, elemI)
            {
                 result[fc[elemI]] += coeffs[elemI]*(pnf[elemI]);
            }
        }
        else
        {
            forAll(fc, elemI)
            {
                result[fc[elemI]] -= coeffs[elemI]*(pnf[elemI]);
            }
        }
    }
}


template<class Type>
void jumpOverlapGgiFvPatchField<Type>::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes,
    const bool switchToLhs
) const
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
