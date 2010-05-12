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

#include "solidWallTemperatureCoupledFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidWallTemperatureCoupledFvPatchScalarField::
solidWallTemperatureCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    coupleManager_(p),
    KName_("undefined-K")
{}


Foam::solidWallTemperatureCoupledFvPatchScalarField::
solidWallTemperatureCoupledFvPatchScalarField
(
    const solidWallTemperatureCoupledFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    coupleManager_(ptf.coupleManager_),
    KName_(ptf.KName_)
{}


Foam::solidWallTemperatureCoupledFvPatchScalarField::
solidWallTemperatureCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    coupleManager_(p, dict),
    KName_(dict.lookup("K"))
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        evaluate();
    }
}


Foam::solidWallTemperatureCoupledFvPatchScalarField::
solidWallTemperatureCoupledFvPatchScalarField
(
    const solidWallTemperatureCoupledFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wtcsf, iF),
    coupleManager_(wtcsf.coupleManager_),
    KName_(wtcsf.KName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidWallTemperatureCoupledFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& neighbourField =
        coupleManager_.neighbourPatchField<scalar>();

    operator==(neighbourField);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::solidWallTemperatureCoupledFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    coupleManager_.writeEntries(os);
    os.writeKeyword("K") << KName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


Foam::tmp<Foam::scalarField>
Foam::solidWallTemperatureCoupledFvPatchScalarField::flux() const
{
    const fvPatchScalarField& Kw =
        patch().lookupPatchField<volScalarField, scalar>(KName_);

    const fvPatchScalarField& Tw = *this;

    return Tw.snGrad()*Kw;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    solidWallTemperatureCoupledFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
