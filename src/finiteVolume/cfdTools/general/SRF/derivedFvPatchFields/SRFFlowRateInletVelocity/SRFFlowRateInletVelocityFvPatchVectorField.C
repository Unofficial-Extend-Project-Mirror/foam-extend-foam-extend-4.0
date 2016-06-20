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

#include "SRFFlowRateInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "SRFModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SRFFlowRateInletVelocityFvPatchVectorField::
SRFFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    flowRate_(0),
    phiName_("phi"),
    rhoName_("rho"),
    relative_(0)
{}


Foam::SRFFlowRateInletVelocityFvPatchVectorField::
SRFFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    flowRate_(readScalar(dict.lookup("flowRate"))),
    phiName_("phi"),
    rhoName_("rho"),
    relative_(dict.lookup("relative"))
{
    if (dict.found("phi"))
    {
        dict.lookup("phi") >> phiName_;
    }

    if (dict.found("rho"))
    {
        dict.lookup("rho") >> rhoName_;
    }
}


Foam::SRFFlowRateInletVelocityFvPatchVectorField::
SRFFlowRateInletVelocityFvPatchVectorField
(
    const SRFFlowRateInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    flowRate_(ptf.flowRate_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    relative_(ptf.relative_)
{}


Foam::SRFFlowRateInletVelocityFvPatchVectorField::
SRFFlowRateInletVelocityFvPatchVectorField
(
    const SRFFlowRateInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    flowRate_(ptf.flowRate_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    relative_(ptf.relative_)
{}


Foam::SRFFlowRateInletVelocityFvPatchVectorField::
SRFFlowRateInletVelocityFvPatchVectorField
(
    const SRFFlowRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    flowRate_(ptf.flowRate_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    relative_(ptf.relative_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SRFFlowRateInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar avgU = -flowRate_/gSum(patch().magSf());

    vectorField n = patch().nf();

    const surfaceScalarField& phi = db().lookupObject<surfaceScalarField>
    (
        phiName_
    );

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        // Volumetric flow-rate

        // If relative, include the effect of the SRF
        if (relative_)
        {
            // Get reference to the SRF model
            const SRF::SRFModel& srf =
                db().lookupObject<SRF::SRFModel>("SRFProperties");

            // Determine patch velocity due to SRF
            const vectorField SRFSurfaceNormalVelocity =
                srf.velocity(patch().Cf());

            operator==(n*avgU - SRFSurfaceNormalVelocity);
        }
        else
        {
            operator==(n*avgU);
        }
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        // Mass flow-rate

        const fvPatchScalarField& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        // If relative, include the effect of the SRF
        if (relative_)
        {
            // Get reference to the SRF model
            const SRF::SRFModel& srf =
                db().lookupObject<SRF::SRFModel>("SRFProperties");

            // Determine patch velocity due to SRF
            const vectorField SRFSurfaceNormalVelocity =
                srf.velocity(patch().Cf());

            operator==(n*avgU/rhop - SRFSurfaceNormalVelocity);
        }
        else
        {
            operator==(n*avgU/rhop);
        }
    }
    else
    {
        FatalErrorIn
        (
            "SRFFlowRateInletVelocityFvPatchVectorField::updateCoeffs()"
        )   << "dimensions of " << phiName_ << " are incorrect" << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << nl << exit(FatalError);
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::SRFFlowRateInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    os.writeKeyword("flowRate") << flowRate_
        << token::END_STATEMENT << nl;

    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }

    if (rhoName_ != "rho")
    {
        os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    }

    os.writeKeyword("relative") << relative_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       SRFFlowRateInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
