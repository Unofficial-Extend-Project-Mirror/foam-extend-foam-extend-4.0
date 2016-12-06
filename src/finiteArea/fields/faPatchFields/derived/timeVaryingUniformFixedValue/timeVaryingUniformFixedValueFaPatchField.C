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

#include "timeVaryingUniformFixedValueFaPatchField.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::timeVaryingUniformFixedValueFaPatchField<Type>::
timeVaryingUniformFixedValueFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    fixedValueFaPatchField<Type>(p, iF),
    timeSeries_()
{}


template<class Type>
Foam::timeVaryingUniformFixedValueFaPatchField<Type>::
timeVaryingUniformFixedValueFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFaPatchField<Type>(p, iF),
    timeSeries_(dict)
{
   if (dict.found("value"))
   {
       faPatchField<Type>::operator==(Field<Type>("value", dict, p.size()));
   }
   else
   {
       updateCoeffs();
   }
}


template<class Type>
Foam::timeVaryingUniformFixedValueFaPatchField<Type>::
timeVaryingUniformFixedValueFaPatchField
(
    const timeVaryingUniformFixedValueFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    fixedValueFaPatchField<Type>(ptf, p, iF, mapper),
    timeSeries_(ptf.timeSeries_)
{}


template<class Type>
Foam::timeVaryingUniformFixedValueFaPatchField<Type>::
timeVaryingUniformFixedValueFaPatchField
(
    const timeVaryingUniformFixedValueFaPatchField<Type>& ptf
)
:
    fixedValueFaPatchField<Type>(ptf),
    timeSeries_(ptf.timeSeries_)
{}


template<class Type>
Foam::timeVaryingUniformFixedValueFaPatchField<Type>::
timeVaryingUniformFixedValueFaPatchField
(
    const timeVaryingUniformFixedValueFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    fixedValueFaPatchField<Type>(ptf, iF),
    timeSeries_(ptf.timeSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::timeVaryingUniformFixedValueFaPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    faPatchField<Type>::operator==
    (
        timeSeries_(this->db().time().timeOutputValue())
    );
    fixedValueFaPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::timeVaryingUniformFixedValueFaPatchField<Type>::write
(
    Ostream& os
) const
{
    faPatchField<Type>::write(os);
    timeSeries_.write(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
