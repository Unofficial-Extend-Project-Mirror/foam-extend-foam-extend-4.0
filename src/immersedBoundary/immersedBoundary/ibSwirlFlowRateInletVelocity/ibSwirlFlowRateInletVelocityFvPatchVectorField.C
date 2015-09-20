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

#include "ibSwirlFlowRateInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
ibSwirlFlowRateInletVelocityFvPatchVectorField::
ibSwirlFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    flowRate_(0),
    phiName_("phi"),
    rhoName_("rho"),
    rpm_(0)
{
    calcGeom();
}


Foam::
ibSwirlFlowRateInletVelocityFvPatchVectorField::
ibSwirlFlowRateInletVelocityFvPatchVectorField
(
    const ibSwirlFlowRateInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    flowRate_(ptf.flowRate_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    rpm_(ptf.rpm_)
{
    calcGeom();
}


Foam::
ibSwirlFlowRateInletVelocityFvPatchVectorField::
ibSwirlFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    flowRate_(readScalar(dict.lookup("flowRate"))),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    rpm_(readScalar(dict.lookup("rpm")))
{
    calcGeom();
}


Foam::
ibSwirlFlowRateInletVelocityFvPatchVectorField::
ibSwirlFlowRateInletVelocityFvPatchVectorField
(
    const ibSwirlFlowRateInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    flowRate_(ptf.flowRate_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    rpm_(ptf.rpm_)
{
    calcGeom();
}


Foam::
ibSwirlFlowRateInletVelocityFvPatchVectorField::
ibSwirlFlowRateInletVelocityFvPatchVectorField
(
    const ibSwirlFlowRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    flowRate_(ptf.flowRate_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    rpm_(ptf.rpm_)
{
    calcGeom();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ibSwirlFlowRateInletVelocityFvPatchVectorField::calcGeom()
{
    // Lookup cellIbMask
    const fvPatchScalarField& patchCellIbMask =
        lookupPatchField<volScalarField, scalar>("cellIbMask");

    scalarField faceMask = patchCellIbMask.patchInternalField();

    totArea_ = gSum(faceMask*patch().magSf());
    avgCenter_ = gSum(faceMask*patch().Cf()*patch().magSf())/totArea_;
    avgNormal_ = gSum(faceMask*patch().Sf())/totArea_;
}


void Foam::ibSwirlFlowRateInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar avgU = -flowRate_/totArea_;

    // Update angular velocity - convert [rpm] to [rad/s]
    vectorField tangentialVelocity =
        (rpm_*mathematicalConstant::pi/30.0)
      * ((patch().Cf() - avgCenter_) ^ avgNormal_);

    vectorField n = patch().nf();

    vectorField& U = *this;

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    // Lookup cellIbMask
    const fvPatchScalarField& patchCellIbMask =
        lookupPatchField<volScalarField, scalar>("cellIbMask");

    scalarField faceMask = patchCellIbMask.patchInternalField();

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        // Volumetric flow-rate
        U = faceMask*(tangentialVelocity + n*avgU);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        // Mass flow-rate
        U = faceMask*(tangentialVelocity + n*avgU/rhop);
    }
    else
    {
        FatalErrorIn
        (
            "ibSwirlFlowRateInletVelocityFvPatchVectorField::updateCoeffs()"
        )   << "dimensions of " << phiName_ << " are incorrect" << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << nl << exit(FatalError);
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::ibSwirlFlowRateInletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("flowRate") << flowRate_ << token::END_STATEMENT << nl;
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    os.writeKeyword("rpm") << rpm_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       ibSwirlFlowRateInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
