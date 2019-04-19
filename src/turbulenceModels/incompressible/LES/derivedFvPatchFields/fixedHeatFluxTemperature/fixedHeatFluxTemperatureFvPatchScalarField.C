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
    Prt_(0.0),
    rhoRef_(0.0),
    c_(0.0)
{}


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
    Pr_(0.0),     // Initialized in constructor body
    Prt_(0.0),    // Initialized in constructor body
    rhoRef_(0.0), // Initialized in constructor body
    c_(0.0)       // Initialized in constructor body
{
    // Read the gradient entry from the dictionary
    if (dict.found("gradient"))
    {
        this->gradient() = scalarField("gradient", dict, p.size());
    }

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

    // Read in all the necessary data as dimensioned scalars
    const dimensionedScalar PrDim =
        dimensionedScalar(transportProperties.lookup("Pr"));

    const dimensionedScalar PrtDim =
        dimensionedScalar(transportProperties.lookup("Prt"));

    const dimensionedScalar rhoRefDim =
        dimensionedScalar(transportProperties.lookup("rhoRef"));

    const dimensionedScalar cDim =
        dimensionedScalar(transportProperties.lookup("c"));

    // Perform sanity checks for dimensions
    if (PrDim.dimensions() != dimless)
    {
        FatalIOErrorIn
        (
            "fixedHeatFluxTemperatureFvPatchScalarField::"
            "fixedHeatFluxTemperatureFvPatchScalarField"
            "\n("
            "\n    const fvPatch& p,"
            "\n    const DimensionedField<scalar, volMesh>& iF,"
            "\n    const dictionary& dict"
            "\n)",
            dict
        ) << "Wrong dimensions for Prandtl number (Pr) detected: "
          << PrDim.dimensions()
          << nl
          << "They should be: " << dimless
          << abort(FatalIOError);
    }

    if (PrtDim.dimensions() != dimless)
    {
        FatalIOErrorIn
        (
            "fixedHeatFluxTemperatureFvPatchScalarField::"
            "fixedHeatFluxTemperatureFvPatchScalarField"
            "\n("
            "\n    const fvPatch& p,"
            "\n    const DimensionedField<scalar, volMesh>& iF,"
            "\n    const dictionary& dict"
            "\n)",
            dict
        ) << "Wrong dimensions for turbulent Prandtl number (Prt) detected: "
          << PrtDim.dimensions()
          << nl
          << "They should be: " << dimless
          << abort(FatalIOError);
    }

    if (rhoRefDim.dimensions() != dimDensity)
    {
        FatalIOErrorIn
        (
            "fixedHeatFluxTemperatureFvPatchScalarField::"
            "fixedHeatFluxTemperatureFvPatchScalarField"
            "\n("
            "\n    const fvPatch& p,"
            "\n    const DimensionedField<scalar, volMesh>& iF,"
            "\n    const dictionary& dict"
            "\n)",
            dict
        ) << "Wrong dimensions for reference density (rhoRef) detected: "
          << rhoRefDim.dimensions()
          << nl
          << "They should be: " << dimDensity
          << abort(FatalIOError);
    }

    if (cDim.dimensions() != dimSpecificHeatCapacity)
    {
        FatalIOErrorIn
        (
            "fixedHeatFluxTemperatureFvPatchScalarField::"
            "fixedHeatFluxTemperatureFvPatchScalarField"
            "\n("
            "\n    const fvPatch& p,"
            "\n    const DimensionedField<scalar, volMesh>& iF,"
            "\n    const dictionary& dict"
            "\n)",
            dict
        ) << "Wrong dimensions for specific heat capacity (c) detected: "
          << cDim.dimensions()
          << nl
          << "They should be: " << dimSpecificHeatCapacity
          << abort(FatalIOError);
    }

    // Store values in data members
    Pr_ = PrDim.value();
    Prt_ = PrtDim.value();
    rhoRef_ = rhoRefDim.value();
    c_ = cDim.value();
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
    Prt_(ptf.Prt_),
    rhoRef_(ptf.rhoRef_),
    c_(ptf.c_)
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
    Prt_(ptf.Prt_),
    rhoRef_(ptf.rhoRef_),
    c_(ptf.c_)
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
    Prt_(ptf.Prt_),
    rhoRef_(ptf.rhoRef_),
    c_(ptf.c_)
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

    if (mesh.foundObject<incompressible::turbulenceModel>("turbulenceModel"))
    {
        const incompressible::turbulenceModel& turb =
            mesh.lookupObject
            <
                incompressible::turbulenceModel
            >("turbulenceModel");

        // Calculate effective kappa at the patch
        const scalarField kappaEffp =
            turb.nu().boundaryField()[patchID]/Pr_
          + turb.nut()().boundaryField()[patchID]/Prt_;

        // Calculate gradient at the boundary
        this->gradient() = heatFlux_/(kappaEffp*rhoRef_*c_);
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
