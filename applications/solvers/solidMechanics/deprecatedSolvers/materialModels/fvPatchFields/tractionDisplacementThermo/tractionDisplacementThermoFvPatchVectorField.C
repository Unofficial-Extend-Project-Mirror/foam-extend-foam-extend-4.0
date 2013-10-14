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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "tractionDisplacementThermoFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "rheologyModel.H"
#include "thermalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tractionDisplacementThermoFvPatchVectorField::
tractionDisplacementThermoFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    UName_("undefined"),
    TName_("undefined"),
    rheologyName_("undefined"),
    thermoName_("undefined"),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


tractionDisplacementThermoFvPatchVectorField::
tractionDisplacementThermoFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    UName_(dict.lookup("U")),
    TName_(dict.lookup("T")),
    rheologyName_(dict.lookup("rheology")),
    thermoName_(dict.lookup("thermo")),
    traction_("traction", dict, p.size()),
    pressure_("pressure", dict, p.size())
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


tractionDisplacementThermoFvPatchVectorField::
tractionDisplacementThermoFvPatchVectorField
(
    const tractionDisplacementThermoFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    UName_(tdpvf.UName_),
    TName_(tdpvf.TName_),
    rheologyName_(tdpvf.rheologyName_),
    thermoName_(tdpvf.thermoName_),
    traction_(tdpvf.traction_, mapper),
    pressure_(tdpvf.pressure_, mapper)
{}


tractionDisplacementThermoFvPatchVectorField::
tractionDisplacementThermoFvPatchVectorField
(
    const tractionDisplacementThermoFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    UName_(tdpvf.UName_),
    TName_(tdpvf.TName_),
    rheologyName_(tdpvf.rheologyName_),
    thermoName_(tdpvf.thermoName_),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


tractionDisplacementThermoFvPatchVectorField::
tractionDisplacementThermoFvPatchVectorField
(
    const tractionDisplacementThermoFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    UName_(tdpvf.UName_),
    TName_(tdpvf.TName_),
    rheologyName_(tdpvf.rheologyName_),
    thermoName_(tdpvf.thermoName_),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionDisplacementThermoFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void tractionDisplacementThermoFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const tractionDisplacementThermoFvPatchVectorField& dmptf =
        refCast<const tractionDisplacementThermoFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void tractionDisplacementThermoFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Looking up rheology
    const rheologyModel& rheology =
        this->db().objectRegistry::
        lookupObject<rheologyModel>(rheologyName_);

    const scalarField mu = rheology.mu()().boundaryField()[patch().index()];
    const scalarField lambda =
        rheology.lambda()().boundaryField()[patch().index()];

    vectorField n = patch().nf();

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>("grad(" +UName_ + ")");

    // Thermal component

    // Looking up thermo
    const thermalModel& thermo =
        this->db().objectRegistry::lookupObject<thermalModel>(thermoName_);

    const fvPatchField<scalar>& T =
        patch().lookupPatchField<volScalarField, scalar>(TName_);

    const scalarField rhoThreeKalpha =
        rheology.rho()().boundaryField()[patch().index()]*
        rheology.threeK()().boundaryField()[patch().index()]*
        thermo.alpha()().boundaryField()[patch().index()];

    const scalarField T0 = thermo.T0()().boundaryField()[patch().index()];

    gradient() =
    (
        (traction_ - (pressure_)*n)
      - (n & (mu*gradU.T() - (mu + lambda)*gradU))
      - n*lambda*tr(gradU)
      + n*rhoThreeKalpha*(T - T0)
    )/(2.0*mu + lambda);

    fixedGradientFvPatchVectorField::updateCoeffs();
}


// Write
void tractionDisplacementThermoFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("U") << UName_ << token::END_STATEMENT << nl;
    os.writeKeyword("T") << TName_ << token::END_STATEMENT << nl;
    os.writeKeyword("rheology") << rheologyName_ << token::END_STATEMENT << nl;
    os.writeKeyword("thermo") << thermoName_ << token::END_STATEMENT << nl;
    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    tractionDisplacementThermoFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
