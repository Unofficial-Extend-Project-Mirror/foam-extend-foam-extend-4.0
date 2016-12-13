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

Author
    Ilaria De Dominicis, General Electric Power, (March 2016)

Contributor
    Hrvoje Jasak, Wikki Ltd.

GE CONFIDENTIAL INFORMATION 2016 General Electric Company. All Rights Reserved

\*---------------------------------------------------------------------------*/

#include "mixingPlaneEnthalpyJumpFvPatchField.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
mixingPlaneEnthalpyJumpFvPatchField<Type>::mixingPlaneEnthalpyJumpFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    jumpMixingPlaneFvPatchField<Type>(p, iF),
    rotating_(false),
    jump_(this->size(), pTraits<Type>::zero)
{}


template<class Type>
mixingPlaneEnthalpyJumpFvPatchField<Type>::mixingPlaneEnthalpyJumpFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    jumpMixingPlaneFvPatchField<Type>(p, iF),
    rotating_(dict.lookup("rotating")),
    jump_(this->size(), pTraits<Type>::zero)
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        this->evaluate(Pstream::blocking);
    }
}


template<class Type>
mixingPlaneEnthalpyJumpFvPatchField<Type>::mixingPlaneEnthalpyJumpFvPatchField
(
    const mixingPlaneEnthalpyJumpFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    jumpMixingPlaneFvPatchField<Type>(ptf, p, iF, mapper),
    rotating_(ptf.rotating_),
    jump_(ptf.jump_, mapper)
{}


template<class Type>
mixingPlaneEnthalpyJumpFvPatchField<Type>::mixingPlaneEnthalpyJumpFvPatchField
(
    const mixingPlaneEnthalpyJumpFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    jumpMixingPlaneFvPatchField<Type>(ptf, iF),
    rotating_(ptf.rotating_),
    jump_(ptf.jump_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void mixingPlaneEnthalpyJumpFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    jumpMixingPlaneFvPatchField<Type>::autoMap(m);
    jump_.autoMap(m);
}


template<class Type>
void mixingPlaneEnthalpyJumpFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    jumpMixingPlaneFvPatchField<Type>::rmap(ptf, addr);

    // rmap jump
    const mixingPlaneEnthalpyJumpFvPatchField<Type>& ejPtf =
        refCast<const mixingPlaneEnthalpyJumpFvPatchField<Type> >(ptf);

    jump_.rmap(ejPtf.jump_, addr);
}


template<class Type>
void mixingPlaneEnthalpyJumpFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("patchType")
        << mixingPlaneFvPatch::typeName << token::END_STATEMENT << nl;
    os.writeKeyword("rotating")
        << rotating_ << token::END_STATEMENT << nl;

    IOstream::streamFormat fmt0 = os.format(IOstream::ASCII);
    os.format(fmt0);

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
