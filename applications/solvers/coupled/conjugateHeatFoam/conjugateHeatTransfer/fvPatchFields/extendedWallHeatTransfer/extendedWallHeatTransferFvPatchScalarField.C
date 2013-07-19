/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "extendedWallHeatTransferFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "radiationConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedWallHeatTransferFvPatchScalarField::extendedWallHeatTransferFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    Tinf_(0.0),
    hc_(0.0),
    alpha_(0.0),
    KName_("undefined-K"),
    radiation_(false)
{}


Foam::extendedWallHeatTransferFvPatchScalarField::extendedWallHeatTransferFvPatchScalarField
(
    const extendedWallHeatTransferFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(p, iF),
    Tinf_(ptf.Tinf_),
    hc_(ptf.hc_),
    alpha_(ptf.alpha_),
    radSources_(ptf.radSources_),
    KName_(ptf.KName_),
    radiation_(ptf.radiation_)
{}


Foam::extendedWallHeatTransferFvPatchScalarField::extendedWallHeatTransferFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    Tinf_(readScalar(dict.lookup("Tinf"))),
    hc_(readScalar(dict.lookup("hc"))),
    alpha_(readScalar(dict.lookup("alpha"))),
    KName_(dict.lookup("K")),
    radiation_(readBool(dict.lookup("radiation")))
{
    refValue() = Tinf_;
    refGrad() = 0.0;
    valueFraction() = 0.0;

    // Construct radiation sources
    PtrList<entry> entries(dict.lookup("radiationSources"));
    radSources_.setSize(entries.size());

    forAll(entries, entryI)
    {
        radSources_.set
        (
            entryI,
            externalRadiationSource::New
            (
                entries[entryI].keyword(),
                entries[entryI].dict(),
                p
            )
        );
    }

    fvPatchField<scalar>::operator=
    (
        scalarField("value", dict, p.size())
    );
}


Foam::extendedWallHeatTransferFvPatchScalarField::extendedWallHeatTransferFvPatchScalarField
(
    const extendedWallHeatTransferFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    Tinf_(tppsf.Tinf_),
    hc_(tppsf.hc_),
    alpha_(tppsf.alpha_),
    radSources_(tppsf.radSources_),
    KName_(tppsf.KName_),
    radiation_(tppsf.radiation_)
{}


Foam::extendedWallHeatTransferFvPatchScalarField::extendedWallHeatTransferFvPatchScalarField
(
    const extendedWallHeatTransferFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    Tinf_(tppsf.Tinf_),
    hc_(tppsf.hc_),
    alpha_(tppsf.alpha_),
    radSources_(tppsf.radSources_),
    KName_(tppsf.KName_),
    radiation_(tppsf.radiation_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::extendedWallHeatTransferFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatch& p = patch();

    scalarField Tc = patchInternalField();
    const scalarField& Tb = *this;

    const scalarField& kappaEff =
        p.lookupPatchField<volScalarField, scalar>(KName_);

    scalarField term = kappaEff*patch().deltaCoeffs();

    scalarField Qri(p.size(), 0);
    if (radiation())
    {
        Qri = p.lookupPatchField<volScalarField, scalar>("Qr");
    }

    scalarField Qrio = 4.0*radiation::sigmaSB.value()*pow4(Tb);

    scalarField Two = Tb - (term*(patchInternalField() - Tb) + Qri)/hc_;

    scalarField Qro(p.size(), 0.0);
    forAll(radSources_, rsI)
    {
        Qro += radSources_[rsI].q(Two);
    }

    scalarField term2 = Tb*hc_*alpha_ + hc_*Qrio + alpha_*Qrio;

    //Info << "Tc = " << Tc << endl;
    //Info << "Tb = "<< (scalarField) *this << endl;
    //Info << "Two = "<< Two << endl;
    //Info << "Qri = "<< Qri << endl;
    //Info << "q = "<< q << endl;

    //Info << "q1 = " << hc_*(*this - Two) << endl;
    //Info << "q2 = " << Qro + alpha_*(Two - Tinf_) << endl;
    //Info << "q3 = "<< term*(Tc - *this) + Qri << endl;

    valueFraction() = 1.0 / (1.0 + Tb*term*(hc_ + alpha_)/term2);
    refValue() =
        Tb*(hc_*alpha_*Tinf_ - hc_*Qro + (hc_ + alpha_)*(Qri + Qrio))/term2;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::extendedWallHeatTransferFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("K") << KName_ << token::END_STATEMENT << nl;
    os.writeKeyword("radiation") << radiation_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tinf") << Tinf_ << token::END_STATEMENT << nl;
    os.writeKeyword("hc") << hc_ << token::END_STATEMENT << nl;
    os.writeKeyword("alpha") << alpha_ << token::END_STATEMENT << nl;

    os  << nl << indent << "radiationSources" << nl
        << indent << token::BEGIN_LIST << incrIndent << nl;
    
    forAll(radSources_, rsI)
    {
        os  << indent << radSources_[rsI].name() << nl
            << indent << token::BEGIN_BLOCK << nl << incrIndent;
        radSources_[rsI].write(os);
        os  << decrIndent << indent << token::END_BLOCK << endl;
    }

    os  << decrIndent << indent
        << token::END_LIST << token::END_STATEMENT << nl << nl;

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        extendedWallHeatTransferFvPatchScalarField
    );
}

// ************************************************************************* //
