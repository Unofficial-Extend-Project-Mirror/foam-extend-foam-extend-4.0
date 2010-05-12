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

#include "engineTimeVaryingTotalPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "interpolateXY.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

engineTimeVaryingTotalPressureFvPatchScalarField::engineTimeVaryingTotalPressureFvPatchScalarField
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
    p0_(p.size(), 0.0),
    timeDataPtr_(NULL),
    engineDB_((refCast<const engineTime>(this->db().time())))
{}


engineTimeVaryingTotalPressureFvPatchScalarField::engineTimeVaryingTotalPressureFvPatchScalarField
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
    p0_("p0", dict, p.size()),
    timeDataFileName_(fileName(dict.lookup("timeDataFileName")).expand()),
    timeDataPtr_(NULL),
    engineDB_((refCast<const engineTime>(this->db().time())))
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
        updateCoeffs();
        fvPatchField<scalar>::operator=(p0_);
    }
}


engineTimeVaryingTotalPressureFvPatchScalarField::engineTimeVaryingTotalPressureFvPatchScalarField
(
    const engineTimeVaryingTotalPressureFvPatchScalarField& ptf,
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
    p0_(ptf.p0_, mapper),
    timeDataFileName_(ptf.timeDataFileName_),
    timeDataPtr_(NULL),
    engineDB_((refCast<const engineTime>(this->db().time())))
{}


engineTimeVaryingTotalPressureFvPatchScalarField::engineTimeVaryingTotalPressureFvPatchScalarField
(
    const engineTimeVaryingTotalPressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    p0_(tppsf.p0_),
    timeDataFileName_(tppsf.timeDataFileName_),
    timeDataPtr_(NULL),
    engineDB_((refCast<const engineTime>(this->db().time())))
{}


engineTimeVaryingTotalPressureFvPatchScalarField::engineTimeVaryingTotalPressureFvPatchScalarField
(
    const engineTimeVaryingTotalPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    p0_(tppsf.p0_),
    timeDataFileName_(tppsf.timeDataFileName_),
    timeDataPtr_(NULL),
    engineDB_((refCast<const engineTime>(this->db().time())))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void engineTimeVaryingTotalPressureFvPatchScalarField::checkTable()
{
    if (!timeDataPtr_.valid())
    {
        timeDataPtr_.reset
        (
            new graph("title", "x", "y", IFstream(timeDataFileName_)())
        );
    }

//    if (this->db().time().value() < min(timeDataPtr_().x()))
    if (engineDB_.theta() < min(timeDataPtr_().x()))
    {
        WarningIn
        (
            "engineTimeVaryingTotalPressureFvPatchScalarField::updateCoeffs()"
        )   << "current time (" << engineDB_.theta()
            << ") is less than the minimum in the data table ("
            << min(timeDataPtr_().x()) << ')' << endl
            << "    Continuing with the value for the smallest time"
            << endl;
    }

//    if (this->db().time().value() > max(timeDataPtr_().x()))
    if (engineDB_.theta() < min(timeDataPtr_().x()))
    {
        WarningIn
        (
            "engineTimeVaryingTotalPressureFvPatchScalarField<Type>::updateCoeffs()"
        )   << "current time (" << engineDB_.theta()
            << ") is greater than the maximum in the data table ("
            << max(timeDataPtr_().x()) << ')' << endl
            << "    Continuing with the value for the largest time"
            << endl;
    }
}



void engineTimeVaryingTotalPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
    p0_.autoMap(m);
}


void engineTimeVaryingTotalPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const engineTimeVaryingTotalPressureFvPatchScalarField& tiptf =
        refCast<const engineTimeVaryingTotalPressureFvPatchScalarField>(ptf);

    p0_.rmap(tiptf.p0_, addr);
}


void engineTimeVaryingTotalPressureFvPatchScalarField::updateCoeffs(const vectorField& Up)
{
    if (updated())
    {
        return;
    }
    
    checkTable();

    p0_=
    (
        interpolateXY
        (
            engineDB_.theta(),
            timeDataPtr_().x(),
            timeDataPtr_().y()
        )
    );

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
            "engineTimeVaryingTotalPressureFvPatchScalarField::updateCoeffs()"
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


void engineTimeVaryingTotalPressureFvPatchScalarField::updateCoeffs()
{
    updateCoeffs(patch().lookupPatchField<volVectorField, vector>(UName_));
}


void engineTimeVaryingTotalPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("U") << UName_ << token::END_STATEMENT << nl;
    os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    os.writeKeyword("psi") << psiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("gamma") << gamma_ << token::END_STATEMENT << endl;
    os.writeKeyword("timeDataFileName")
        << timeDataFileName_ << token::END_STATEMENT << nl;
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, engineTimeVaryingTotalPressureFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
