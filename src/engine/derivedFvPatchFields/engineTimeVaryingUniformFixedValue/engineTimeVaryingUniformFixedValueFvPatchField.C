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

#include "engineTimeVaryingUniformFixedValueFvPatchField.H"
#include "graph.H"
#include "IFstream.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
engineTimeVaryingUniformFixedValueFvPatchField<Type>::
engineTimeVaryingUniformFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    timeDataPtr_(NULL),
    engineDB_((refCast<const engineTime>(this->db().time())))
{}


template<class Type>
engineTimeVaryingUniformFixedValueFvPatchField<Type>::
engineTimeVaryingUniformFixedValueFvPatchField
(
    const engineTimeVaryingUniformFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    timeDataFileName_(ptf.timeDataFileName_),
    timeDataPtr_(NULL),
    engineDB_((refCast<const engineTime>(this->db().time())))
{}


template<class Type>
engineTimeVaryingUniformFixedValueFvPatchField<Type>::
engineTimeVaryingUniformFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    timeDataFileName_(fileName(dict.lookup("timeDataFileName")).expand()),
    timeDataPtr_(NULL),
    engineDB_((refCast<const engineTime>(this->db().time())))
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator==(Field<Type>("value", dict, p.size()));
    }
    else
    {
        updateCoeffs();
    }
}


template<class Type>
engineTimeVaryingUniformFixedValueFvPatchField<Type>::
engineTimeVaryingUniformFixedValueFvPatchField
(
    const engineTimeVaryingUniformFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    timeDataFileName_(ptf.timeDataFileName_),
    timeDataPtr_(NULL),
    engineDB_((refCast<const engineTime>(this->db().time())))
{}


template<class Type>
engineTimeVaryingUniformFixedValueFvPatchField<Type>::
engineTimeVaryingUniformFixedValueFvPatchField
(
    const engineTimeVaryingUniformFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    timeDataFileName_(ptf.timeDataFileName_),
    timeDataPtr_(NULL),
    engineDB_((refCast<const engineTime>(this->db().time())))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void engineTimeVaryingUniformFixedValueFvPatchField<Type>::checkTable()
{
    if (!timeDataPtr_.valid())
    {
        timeDataPtr_.reset
        (
            new graph("title", "x", "y", IFstream(timeDataFileName_)())
        );
    }

//    if (this->db().time().value() < min(timeDataPtr_().x()))
    if (engineDB_.theta() < min(timeDataPtr_().x()))
    {
        WarningIn
        (
            "engineTimeVaryingUniformFixedValueFvPatchField<Type>::updateCoeffs()"
        )   << "current time (" << engineDB_.theta()
            << ") is less than the minimum in the data table ("
            << min(timeDataPtr_().x()) << ')' << endl
            << "    Continuing with the value for the smallest time"
            << endl;
    }

//    if (this->db().time().value() > max(timeDataPtr_().x()))
    if (engineDB_.theta() < min(timeDataPtr_().x()))
    {
        WarningIn
        (
            "engineTimeVaryingUniformFixedValueFvPatchField<Type>::updateCoeffs()"
        )   << "current time (" << engineDB_.theta()
            << ") is greater than the maximum in the data table ("
            << max(timeDataPtr_().x()) << ')' << endl
            << "    Continuing with the value for the largest time"
            << endl;
    }
}


template<class Type>
void engineTimeVaryingUniformFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("timeDataFileName")
        << timeDataFileName_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
