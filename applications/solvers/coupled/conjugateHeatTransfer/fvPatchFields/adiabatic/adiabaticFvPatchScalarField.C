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

#include "adiabaticFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "radiationConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adiabaticFvPatchScalarField::
adiabaticFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    KName_("undefined-K")
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::adiabaticFvPatchScalarField::
adiabaticFvPatchScalarField
(
    const adiabaticFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    KName_(ptf.KName_)
{}


Foam::adiabaticFvPatchScalarField::
adiabaticFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    KName_(dict.lookup("K"))
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
}


Foam::adiabaticFvPatchScalarField::
adiabaticFvPatchScalarField
(
    const adiabaticFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    KName_(tppsf.KName_)
{}


Foam::adiabaticFvPatchScalarField::
adiabaticFvPatchScalarField
(
    const adiabaticFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    KName_(tppsf.KName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adiabaticFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& Tp = *this;

    const scalarField& Kw =
        patch().lookupPatchField<volScalarField, scalar>(KName_);

    const scalarField& Qr =
        patch().lookupPatchField<volScalarField, scalar>("Qr");

    scalarField fourQro = 4.0*radiation::sigmaSB.value()*pow4(Tp);

    refGrad() = (Qr + fourQro)/Kw;
    valueFraction() = fourQro/(fourQro + Kw*patch().deltaCoeffs()*Tp);

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::adiabaticFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("K") << KName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        adiabaticFvPatchScalarField
    );
}

// ************************************************************************* //
