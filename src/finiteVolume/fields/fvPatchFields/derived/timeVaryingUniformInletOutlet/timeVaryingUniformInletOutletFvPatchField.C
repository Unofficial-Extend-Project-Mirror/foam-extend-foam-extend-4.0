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

\*---------------------------------------------------------------------------*/

#include "timeVaryingUniformInletOutletFvPatchField.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::timeVaryingUniformInletOutletFvPatchField<Type>::
timeVaryingUniformInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    inletOutletFvPatchField<Type>(p, iF),
    timeSeries_()
{}


template<class Type>
Foam::timeVaryingUniformInletOutletFvPatchField<Type>::
timeVaryingUniformInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchField<Type>(p, iF),
    timeSeries_(dict)
{
    this->refValue() = timeSeries_(this->db().time().timeOutputValue());
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator==(Field<Type>("value", dict, p.size()));
    }
    else
    {
        fvPatchField<Type>::operator=(this->refValue());
    }
}


template<class Type>
Foam::timeVaryingUniformInletOutletFvPatchField<Type>::
timeVaryingUniformInletOutletFvPatchField
(
    const timeVaryingUniformInletOutletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchField<Type>(ptf, p, iF, mapper),
    timeSeries_(ptf.timeSeries_)
{}


template<class Type>
Foam::timeVaryingUniformInletOutletFvPatchField<Type>::
timeVaryingUniformInletOutletFvPatchField
(
    const timeVaryingUniformInletOutletFvPatchField<Type>& ptf
)
:
    inletOutletFvPatchField<Type>(ptf),
    timeSeries_(ptf.timeSeries_)
{}


template<class Type>
Foam::timeVaryingUniformInletOutletFvPatchField<Type>::
timeVaryingUniformInletOutletFvPatchField
(
    const timeVaryingUniformInletOutletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    inletOutletFvPatchField<Type>(ptf, iF),
    timeSeries_(ptf.timeSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::timeVaryingUniformInletOutletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    this->refValue() = timeSeries_(this->db().time().timeOutputValue());
    inletOutletFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::timeVaryingUniformInletOutletFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fvPatchField<Type>::write(os);
    timeSeries_.write(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
