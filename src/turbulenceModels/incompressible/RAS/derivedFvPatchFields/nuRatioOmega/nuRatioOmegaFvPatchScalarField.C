/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "nuRatioOmegaFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nuRatioOmegaFvPatchScalarField::nuRatioOmegaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    nuRatio_(1),
    kName_("undefined-k"),
    phiName_("undefined-phi")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::nuRatioOmegaFvPatchScalarField::nuRatioOmegaFvPatchScalarField
(
    const nuRatioOmegaFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    nuRatio_(ptf.nuRatio_),
    kName_(ptf.kName_),
    phiName_(ptf.phiName_)
{}


Foam::nuRatioOmegaFvPatchScalarField::nuRatioOmegaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    nuRatio_(readScalar(dict.lookup("nuRatio"))),
    kName_(dict.lookupOrDefault<word>("kName", "k")),
    phiName_(dict.lookupOrDefault<word>("phiName", "phi"))
{
    if (nuRatio_< SMALL)
    {
        FatalErrorIn
        (
            "nuRatioOmegaFvPatchScalarField::nuRatioOmegaFvPatchScalarField"
            "(const fvPatch& p, const DimensionedField<scalar, volMesh>& iF, "
            "const dictionary& dict)"
        )   << "Invalid eddy viscosity ratio (nuRatio) specified: " << nuRatio_
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::nuRatioOmegaFvPatchScalarField::nuRatioOmegaFvPatchScalarField
(
    const nuRatioOmegaFvPatchScalarField& ptf
)
:
    inletOutletFvPatchScalarField(ptf),
    nuRatio_(ptf.nuRatio_),
    kName_(ptf.kName_),
    phiName_(ptf.phiName_)
{}


Foam::nuRatioOmegaFvPatchScalarField::nuRatioOmegaFvPatchScalarField
(
    const nuRatioOmegaFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF),
    nuRatio_(ptf.nuRatio_),
    kName_(ptf.kName_),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nuRatioOmegaFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get turbulent kinetic energy field for this patch
    const fvPatchScalarField& kp =
        lookupPatchField<volScalarField, vector>(kName_);

    // Get flux field for this patch
    const fvsPatchScalarField& phip =
        lookupPatchField<surfaceScalarField, scalar>(phiName_);

    // Get RASModel
    const incompressible::RASModel& rasModel =
        this->dimensionedInternalField().mesh().lookupObject
        <
            incompressible::RASModel
        >("RASProperties");

    // Get laminar viscosity for this patch
    const fvPatchScalarField& nup =
        rasModel.nu().boundaryField()[this->patch().index()];

    this->refValue() = kp/(nuRatio_*nup);
    this->valueFraction() = neg(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}


void Foam::nuRatioOmegaFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("nuRatio") << nuRatio_ << token::END_STATEMENT << nl;
    os.writeKeyword("kName") << kName_ << token::END_STATEMENT << nl;
    os.writeKeyword("phiName") << phiName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nuRatioOmegaFvPatchScalarField
    );
}

// ************************************************************************* //
