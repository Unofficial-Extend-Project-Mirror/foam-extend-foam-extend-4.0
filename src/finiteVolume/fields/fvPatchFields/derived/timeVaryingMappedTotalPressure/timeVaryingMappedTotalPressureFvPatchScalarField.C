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

#include "timeVaryingMappedTotalPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeVaryingMappedTotalPressureFvPatchScalarField::
timeVaryingMappedTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    timeVaryingMappedFixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi"),
    rhoName_("none"),
    rho_(1),
    psiName_("none"),
    gamma_(0.0)
{}


Foam::timeVaryingMappedTotalPressureFvPatchScalarField::
timeVaryingMappedTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    timeVaryingMappedFixedValueFvPatchScalarField(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("UName", "U")),
    phiName_(dict.lookupOrDefault<word>("phiName", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rhoName", "none")),
    rho_
    (
        iF.dimensions() == dimPressure/dimDensity && rhoName_ == "rho"
        ? readScalar(dict.lookup("rho"))
        : 1
    ),
    psiName_(dict.lookupOrDefault<word>("psiName", "none")),
    gamma_
    (
        iF.dimensions() == dimPressure && psiName_ != "none"
        ? readScalar(dict.lookup("gamma"))
        : 1
    )
{}


Foam::timeVaryingMappedTotalPressureFvPatchScalarField::
timeVaryingMappedTotalPressureFvPatchScalarField
(
    const timeVaryingMappedTotalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    timeVaryingMappedFixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    rho_(ptf.rho_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_)
{}


Foam::timeVaryingMappedTotalPressureFvPatchScalarField::
timeVaryingMappedTotalPressureFvPatchScalarField
(
    const timeVaryingMappedTotalPressureFvPatchScalarField& tppsf
)
:
    timeVaryingMappedFixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    rho_(tppsf.rho_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_)
{}


Foam::timeVaryingMappedTotalPressureFvPatchScalarField::
timeVaryingMappedTotalPressureFvPatchScalarField
(
    const timeVaryingMappedTotalPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    timeVaryingMappedFixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    rho_(tppsf.rho_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeVaryingMappedTotalPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    timeVaryingMappedFixedValueFvPatchScalarField::autoMap(m);
}


void Foam::timeVaryingMappedTotalPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    timeVaryingMappedFixedValueFvPatchScalarField::rmap(ptf, addr);
}


void Foam::timeVaryingMappedTotalPressureFvPatchScalarField::updateCoeffs
(
    const vectorField& Up
)
{
    if (updated())
    {
        return;
    }

    // Map the total-pressure field onto this field
    timeVaryingMappedFixedValueFvPatchScalarField::updateCoeffs();
    const scalarField p0 = *this;

    const fvsPatchField<scalar>& phip =
        lookupPatchField<surfaceScalarField, scalar>(phiName_);

    if (dimensionedInternalField().dimensions() == dimPressure/dimDensity)
    {
        if (rhoName_ == "none")
        {
            operator==(p0 - 0.5*(1.0 - pos(phip))*magSqr(Up));
        }
        else if (rhoName_ == "rho")
        {
            operator==(p0/rho_ - 0.5*(1.0 - pos(phip))*magSqr(Up));
        }
        else
        {
            FatalErrorIn
            (
                "timeVaryingMappedTotalPressureFvPatchScalarField::"
                "updateCoeffs()"
            )   << " rhoName set inconsistently, rhoName = " << rhoName_
                << ".\n"
                << "    Set rhoName to 'rho' or 'none' depending on the "
                   "definition of total pressure." << nl
                << "    on patch " << this->patch().name()
                << " of field " << this->dimensionedInternalField().name()
                << " in file " << this->dimensionedInternalField().objectPath()
                << exit(FatalError);
        }
    }
    else if (dimensionedInternalField().dimensions() == dimPressure)
    {
        if (rhoName_ != "none")
        {
            const fvPatchField<scalar>& rho =
                lookupPatchField<volScalarField, scalar>(rhoName_);

            operator==(p0 - 0.5*rho*(1.0 - pos(phip))*magSqr(Up));
        }
        else if (psiName_ != "none")
        {
            const fvPatchField<scalar>& psip =
                lookupPatchField<volScalarField, scalar>(psiName_);

            if (gamma_ > 1.0)
            {
                scalar gM1ByG = (gamma_ - 1.0)/gamma_;

                operator==
                (
                    p0
                    /pow
                    (
                        (1.0 + 0.5*psip*gM1ByG*(1.0 - pos(phip))*magSqr(Up)),
                        1.0/gM1ByG
                    )
                );
            }
            else
            {
                operator==(p0/(1.0 + 0.5*psip*(1.0 - pos(phip))*magSqr(Up)));
            }
        }
        else
        {
            FatalErrorIn
            (
                "timeVaryingMappedTotalPressureFvPatchScalarField::"
                "updateCoeffs()"
            )   << " rhoName or psiName set inconsistently, rhoName = "
                << rhoName_ << ", psiName = " << psiName_ << ".\n"
                << "    Set either rhoName or psiName depending on the "
                   "definition of total pressure." << nl
                << "    Set the unused variable(s) to 'none'.\n"
                << "    on patch " << this->patch().name()
                << " of field " << this->dimensionedInternalField().name()
                << " in file " << this->dimensionedInternalField().objectPath()
                << exit(FatalError);
        }
    }
    else
    {
        FatalErrorIn
        (
            "timeVaryingMappedTotalPressureFvPatchScalarField::updateCoeffs()"
        )   << "Incorrect dimensions for pressure "
            << dimensionedInternalField().dimensions()
            << "    Should be either " << dimPressure
            << " for compressible/variable density flow\n"
            << "    or " << dimPressure/dimDensity
            << " for incompressible flow.\n"
            << "    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }
}


void Foam::timeVaryingMappedTotalPressureFvPatchScalarField::updateCoeffs()
{
    updateCoeffs(lookupPatchField<volVectorField, vector>(UName_));
}


void Foam::timeVaryingMappedTotalPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    writeEntryIfDifferent<word>(os, "UName", "U", UName_);
    writeEntryIfDifferent<word>(os, "phiName", "phi", phiName_);
    os.writeKeyword("rhoName") << rhoName_ << token::END_STATEMENT << nl;
    os.writeKeyword("rho") << rho_ << token::END_STATEMENT << nl;
    os.writeKeyword("psiName") << psiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("gamma") << gamma_ << token::END_STATEMENT << nl;
    timeVaryingMappedFixedValueFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        timeVaryingMappedTotalPressureFvPatchScalarField
    );
}

// ************************************************************************* //
