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

#include "fixedFluxPressureFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi"),
    rhoName_("rho"),
    rAUName_("(1|A(" + UName_ + "))"),
    adjoint_(false)
{}


Foam::fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    rAUName_(dict.lookupOrDefault<word>("rAU", "(1|A(" + UName_ + "))")),
    adjoint_(dict.lookupOrDefault<Switch>("adjoint", false))
{
    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fixedFluxPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    rAUName_(ptf.rAUName_),
    adjoint_(ptf.adjoint_)
{}


Foam::fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fixedFluxPressureFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    rAUName_(ptf.rAUName_),
    adjoint_(ptf.adjoint_)
{}


Foam::fixedFluxPressureFvPatchScalarField::fixedFluxPressureFvPatchScalarField
(
    const fixedFluxPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    rAUName_(ptf.rAUName_),
    adjoint_(ptf.adjoint_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedFluxPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<vector>& Up =
        lookupPatchField<volVectorField, vector>(UName_);

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    fvsPatchField<scalar> phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            lookupPatchField<volScalarField, scalar>(rhoName_);

        phip /= rhop;
    }

    // Collect diffusivity from the pressure equation
    scalarField rAp;

    if (this->db().objectRegistry::foundObject<volScalarField>(rAUName_))
    {
        rAp = lookupPatchField<volScalarField, scalar>(rAUName_);
    }
    else if
    (
        this->db().objectRegistry::foundObject<surfaceScalarField>(rAUName_)
    )
    {
        rAp = lookupPatchField<surfaceScalarField, scalar>(rAUName_);
    }
    else
    {
        FatalErrorIn
        (
            "void fixedFluxPressureFvPatchScalarField::updateCoeffs()"
        )   << " Field " << rAUName_ << " not found for patch "
            << patch().name() << " and field "
            << dimensionedInternalField().name()
            << abort(FatalError);
    }

    if (adjoint_)
    {
        gradient() = ((patch().Sf() & Up) - phip)/patch().magSf()/rAp;
    }
    else
    {
        gradient() = (phip - (patch().Sf() & Up))/patch().magSf()/rAp;
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::fixedFluxPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "rAU", "rAU", rAUName_);
    writeEntryIfDifferent<Switch>(os, "adjoint", false, adjoint_);

    gradient().writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedFluxPressureFvPatchScalarField
    );
}

// ************************************************************************* //
