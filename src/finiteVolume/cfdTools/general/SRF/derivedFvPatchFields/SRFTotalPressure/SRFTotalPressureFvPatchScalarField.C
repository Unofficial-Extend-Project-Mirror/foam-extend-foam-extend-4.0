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

#include "SRFTotalPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "SRFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SRFTotalPressureFvPatchScalarField::SRFTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("undefined"),
    phiName_("undefined"),
    rhoName_("undefined"),
    psiName_("undefined"),
    gamma_(0.0),
    relative_(0),
    p0_(p.size(), 0.0)
{}


SRFTotalPressureFvPatchScalarField::SRFTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_(dict.lookup("U")),
    phiName_(dict.lookup("phi")),
    rhoName_(dict.lookup("rho")),
    psiName_(dict.lookup("psi")),
    gamma_(readScalar(dict.lookup("gamma"))),
    relative_(dict.lookup("relative")),
    p0_("p0", dict, p.size())
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
        fvPatchField<scalar>::operator=(p0_);
    }
}


SRFTotalPressureFvPatchScalarField::SRFTotalPressureFvPatchScalarField
(
    const SRFTotalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    relative_(ptf.relative_),
    p0_(ptf.p0_, mapper)
{}


SRFTotalPressureFvPatchScalarField::SRFTotalPressureFvPatchScalarField
(
    const SRFTotalPressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    relative_(tppsf.relative_),
    p0_(tppsf.p0_)
{}


SRFTotalPressureFvPatchScalarField::SRFTotalPressureFvPatchScalarField
(
    const SRFTotalPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    relative_(tppsf.relative_),
    p0_(tppsf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SRFTotalPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    p0_.autoMap(m);
}


void SRFTotalPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const SRFTotalPressureFvPatchScalarField& tiptf =
        refCast<const SRFTotalPressureFvPatchScalarField>(ptf);

    p0_.rmap(tiptf.p0_, addr);
}


void SRFTotalPressureFvPatchScalarField::updateCoeffs(const vectorField& Up)
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    if (psiName_ == "none" && rhoName_ == "none")
    {
        operator==(p0_ - 0.5*(1.0 - pos(phip))*magSqr(Up));
    }
    else if (rhoName_ == "none")
    {
        const fvPatchField<scalar>& psip =
            patch().lookupPatchField<volScalarField, scalar>(psiName_);

        if (gamma_ > 1.0)
        {
            scalar gM1ByG = (gamma_ - 1.0)/gamma_;

            operator==
            (
                p0_
               /pow
                (
                    (1.0 + 0.5*psip*gM1ByG*(1.0 - pos(phip))*magSqr(Up)),
                    1.0/gM1ByG
                )
            );
        }
        else
        {
            operator==(p0_/(1.0 + 0.5*psip*(1.0 - pos(phip))*magSqr(Up)));
        }
    }
    else if (psiName_ == "none")
    {
        const fvPatchField<scalar>& rho =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        operator==(p0_ - 0.5*rho*(1.0 - pos(phip))*magSqr(Up));
    }
    else
    {
        FatalErrorIn
        (
            "SRFTotalPressureFvPatchScalarField::updateCoeffs()"
        )   << " rho or psi set inconsitently, rho = " << rhoName_
            << ", psi = " << psiName_ << '.' << nl
            << "    Set either rho or psi or neither depending on the "
               "definition of total pressure." << nl
            << "    Set the unused variables to 'none'."
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void SRFTotalPressureFvPatchScalarField::updateCoeffs()
{
    const vectorField& U =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    // If relative, include the effect of the SRF
    if (relative_)
    {
        // Get reference to the SRF model
        const SRF::SRFModel& srf =
            db().lookupObject<SRF::SRFModel>("SRFProperties");

        // Determine patch velocity due to SRF
        const vectorField SRFSurfaceNormalVelocity =
            srf.velocity(patch().Cf());

        updateCoeffs(U - SRFSurfaceNormalVelocity);
    }
    else
    {
        updateCoeffs(U);
    }
}


void SRFTotalPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("U") << UName_ << token::END_STATEMENT << nl;
    os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    os.writeKeyword("psi") << psiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("gamma") << gamma_ << token::END_STATEMENT << endl;
    os.writeKeyword("relative") << relative_ << token::END_STATEMENT << nl;
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, SRFTotalPressureFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
