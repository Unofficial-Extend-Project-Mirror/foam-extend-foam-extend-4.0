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

#include "fixedRhoEFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedRhoEFvPatchScalarField::fixedRhoEFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


fixedRhoEFvPatchScalarField::fixedRhoEFvPatchScalarField
(
    const fixedRhoEFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


fixedRhoEFvPatchScalarField::fixedRhoEFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


fixedRhoEFvPatchScalarField::fixedRhoEFvPatchScalarField
(
    const fixedRhoEFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf)
{}


fixedRhoEFvPatchScalarField::fixedRhoEFvPatchScalarField
(
    const fixedRhoEFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedRhoEFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const dictionary& thermodynamicProperties = db().lookupObject<IOdictionary>
    (
        "thermodynamicProperties"
    );

    dimensionedScalar Cv(thermodynamicProperties.lookup("Cv"));

    const fvPatchScalarField& rhop =
        lookupPatchField<volScalarField, scalar>("rho");

    const fvPatchVectorField& rhoUp =
        lookupPatchField<volVectorField, vector>("rhoU");

    const fvPatchScalarField& Tp =
        lookupPatchField<volScalarField, scalar>("T");

    operator==(rhop*(Cv.value()*Tp + 0.5*magSqr(rhoUp/rhop)));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void fixedRhoEFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, fixedRhoEFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
