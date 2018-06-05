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

#include "fixedHeatFluxTemperatureFvPatchScalarField.H"
#include "incompressible/RAS/RASModel/RASModel.H"
#include "incompressible/LES/LESModel/LESModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedHeatFluxTemperatureFvPatchScalarField::
fixedHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    heatFlux_(p.size(), 0.0),
    Pr_(0.0),
    Prt_(0.0)
{
    this->gradient() = 0.0;
}


fixedHeatFluxTemperatureFvPatchScalarField::
fixedHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    heatFlux_("heatFlux", dict, p.size()),
    Pr_(0.0), // Initialized in constructor body
    Prt_(0.0) // Initialized in constructor body
{
    // Set dummy gradient
    this->gradient() = 0.0;

    // Read the value entry from the dictionary
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        FatalIOErrorIn
        (
            "fixedHeatFluxTemperatureFvPatchScalarField::"
            "fixedHeatFluxTemperatureFvPatchScalarField"
            "("
            "const fvPatch& p,"
            "const DimensionedField<scalar, volMesh>& iF,"
            "const dictionary& dict,"
            "const bool valueRequired"
            ")",
            dict
        )   << "Essential entry 'value' missing"
            << exit(FatalIOError);
    }

    // Get mesh
    const fvMesh& mesh = this->dimensionedInternalField().mesh();

    // Create transportProperties dictionary
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false // do not register in the database
        )
    );

    // Read in Prandtl numbers
    Pr_ = dimensionedScalar(transportProperties.lookup("Pr")).value();
    Prt_ = dimensionedScalar(transportProperties.lookup("Prt")).value();
}


fixedHeatFluxTemperatureFvPatchScalarField::
fixedHeatFluxTemperatureFvPatchScalarField
(
    const fixedHeatFluxTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    heatFlux_(ptf.heatFlux_),
    Pr_(ptf.Pr_),
    Prt_(ptf.Prt_)
{}


fixedHeatFluxTemperatureFvPatchScalarField::
fixedHeatFluxTemperatureFvPatchScalarField
(
    const fixedHeatFluxTemperatureFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    heatFlux_(ptf.heatFlux_),
    Pr_(ptf.Pr_),
    Prt_(ptf.Prt_)
{}


fixedHeatFluxTemperatureFvPatchScalarField::
fixedHeatFluxTemperatureFvPatchScalarField
(
    const fixedHeatFluxTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    heatFlux_(ptf.heatFlux_),
    Pr_(ptf.Pr_),
    Prt_(ptf.Prt_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Calculate the gradient depending on the turbulence model
    const fvMesh& mesh = this->dimensionedInternalField().mesh();

    // Get this patch index
    const label patchID = this->patch().index();

    if (mesh.foundObject<incompressible::RASModel>("RASProperties"))
    {
        const incompressible::RASModel& ras =
            mesh.lookupObject<incompressible::RASModel>("RASProperties");

        // Calculate effective kappa at the patch
        const scalarField kappaEffp =
            ras.nu().boundaryField()[patchID]/Pr_
          + ras.nut()().boundaryField()[patchID]/Prt_;

        this->gradient() = heatFlux_/kappaEffp;
    }
    else if (mesh.foundObject<incompressible::LESModel>("LESProperties"))
    {
        const incompressible::LESModel& les =
            mesh.lookupObject<incompressible::LESModel>("LESProperties");

        // Calculate effective kappa at the patch
        const scalarField kappaEffp =
            les.nu().boundaryField()[patchID]/Pr_
          + les.nut()().boundaryField()[patchID]/Prt_;

        this->gradient() = heatFlux_/kappaEffp;
    }
    else
    {
        FatalErrorIn
        (
            "fixedHeatFluxTemperatureFvPatchScalarField::updateCoeffs()"
        )   << " No valid model for effective kappa calculations."
            << abort(FatalError);
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void fixedHeatFluxTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    heatFlux_.writeEntry("heatFlux", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    fixedHeatFluxTemperatureFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
