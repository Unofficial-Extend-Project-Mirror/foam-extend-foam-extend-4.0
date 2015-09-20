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

#include "engineTimeVaryingUniformFixedValueFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makePatchFields(engineTimeVaryingUniformFixedValue);

// Update the coefficients associated with the patch field
template<>
void engineTimeVaryingUniformFixedValueFvPatchField<scalar>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    checkTable();

    this->operator==
    (
        interpolateXY
        (
//            this->db().time().value(),
            engineDB_.theta(),
            timeDataPtr_().x(),
            timeDataPtr_().y()
        )
    );

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


template<>
void engineTimeVaryingUniformFixedValueFvPatchField<vector>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    checkTable();

    this->operator==
    (
       -interpolateXY
        (
//            this->db().time().value(),
            engineDB_.theta(),
            timeDataPtr_().x(),
            timeDataPtr_().y()
        )*patch().nf()
    );

    fixedValueFvPatchField<vector>::updateCoeffs();
}


template<>
void engineTimeVaryingUniformFixedValueFvPatchField<sphericalTensor>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    checkTable();

    this->operator==
    (
        interpolateXY
        (
//            this->db().time().value(),
            engineDB_.theta(),
            timeDataPtr_().x(),
            timeDataPtr_().y()
        )*sphericalTensor::I
    );

    fixedValueFvPatchField<sphericalTensor>::updateCoeffs();
}


template<>
void engineTimeVaryingUniformFixedValueFvPatchField<symmTensor>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    checkTable();

    this->operator==
    (
        interpolateXY
        (
//            this->db().time().value(),
            engineDB_.theta(),
            timeDataPtr_().x(),
            timeDataPtr_().y()
        )*sqr(patch().nf())
    );

    fixedValueFvPatchField<symmTensor>::updateCoeffs();
}


template<>
void engineTimeVaryingUniformFixedValueFvPatchField<symmTensor4thOrder>::updateCoeffs()
{
    notImplemented
    (
        "engineTimeVaryingUniformFixedValueFvPatchField"
        "<symmTensor4thOrder>::updateCoeffs()"
    );
}


template<>
void engineTimeVaryingUniformFixedValueFvPatchField<diagTensor>::updateCoeffs()
{
    notImplemented
    (
        "engineTimeVaryingUniformFixedValueFvPatchField"
        "<diagTensor>::updateCoeffs()"
    );
}


template<>
void engineTimeVaryingUniformFixedValueFvPatchField<tensor>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    checkTable();

    this->operator==
    (
        interpolateXY
        (
//            this->db().time().value(),
            engineDB_.theta(),
            timeDataPtr_().x(),
            timeDataPtr_().y()
        )*(patch().nf()*patch().nf())
    );

    fixedValueFvPatchField<tensor>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
