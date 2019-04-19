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

#include "kNonEqWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "RASModel.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void kNonEqWallFunctionFvPatchScalarField::checkType()
{
    if (!this->patch().isWall())
    {
        FatalErrorIn("kNonEqWallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << this->patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << this->patch().type()
            << nl << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kNonEqWallFunctionFvPatchScalarField::kNonEqWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF

)
:
    zeroGradientFvPatchScalarField(p, iF),
    UName_("U"),
    kName_("k"),
    epsilonName_("epsilon"),
    GName_("RASModel::G"),
    nuName_("nu"),
    nutName_("nut"),
    Cmu_(0.09),
    kappa_(0.41)
{
    checkType();
}


kNonEqWallFunctionFvPatchScalarField::kNonEqWallFunctionFvPatchScalarField
(
    const kNonEqWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper

)
:
    zeroGradientFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    kName_(ptf.kName_),
    epsilonName_(ptf.epsilonName_),
    GName_(ptf.GName_),
    nuName_(ptf.nuName_),
    nutName_(ptf.nutName_),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_)

{
    checkType();
}


kNonEqWallFunctionFvPatchScalarField::kNonEqWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    kName_(dict.lookupOrDefault<word>("k", "k")),
    epsilonName_(dict.lookupOrDefault<word>("epsilon", "epsilon")),
    GName_(dict.lookupOrDefault<word>("G", "RASModel::G")),
    nuName_(dict.lookupOrDefault<word>("nu", "nu")),
    nutName_(dict.lookupOrDefault<word>("nut", "nut")),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41))

{
    checkType();
}


kNonEqWallFunctionFvPatchScalarField::kNonEqWallFunctionFvPatchScalarField
(
    const kNonEqWallFunctionFvPatchScalarField& tkqrwfpf
)
:
    zeroGradientFvPatchScalarField(tkqrwfpf),
    UName_(tkqrwfpf.UName_),
    kName_(tkqrwfpf.kName_),
    epsilonName_(tkqrwfpf.epsilonName_),
    GName_(tkqrwfpf.GName_),
    nuName_(tkqrwfpf.nuName_),
    nutName_(tkqrwfpf.nutName_),
    Cmu_(tkqrwfpf.Cmu_),
    kappa_(tkqrwfpf.kappa_)
{
    checkType();
}


kNonEqWallFunctionFvPatchScalarField::kNonEqWallFunctionFvPatchScalarField
(
    const kNonEqWallFunctionFvPatchScalarField& tkqrwfpf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(tkqrwfpf, iF),
    UName_(tkqrwfpf.UName_),
    kName_(tkqrwfpf.kName_),
    epsilonName_(tkqrwfpf.epsilonName_),
    GName_(tkqrwfpf.GName_),
    nuName_(tkqrwfpf.nuName_),
    nutName_(tkqrwfpf.nutName_),
    Cmu_(tkqrwfpf.Cmu_),
    kappa_(tkqrwfpf.kappa_)
{
    checkType();
}


void kNonEqWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // If G field is not present, execute zero gradient evaluation
    // HJ, 20/Mar/2011
    if ( !db().foundObject<volScalarField>(GName_) )

    {
        InfoIn("void kNonEqWallFunctionFvPatchScalarField::updateCoeffs()")
            << "Cannot access " << GName_ << " field for patch "
            << patch().name() << ".  Evaluating as zeroGradient"
            << endl;

        fvPatchScalarField::updateCoeffs();

        zeroGradientFvPatchScalarField::evaluate();


        return;
    }

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalarField& y = rasModel.y()[patch().index()];

    const scalar Cmu25 = pow(Cmu_, 0.25);
    const scalar Cmu75 = pow(Cmu_, 0.75);

    volScalarField& G = const_cast<volScalarField&>
        (db().lookupObject<volScalarField>(GName_));
    scalarField& GIn = G.internalField();

    volScalarField& epsilon = const_cast<volScalarField&>
        (db().lookupObject<volScalarField>(epsilonName_));
    scalarField& epsilonIn = epsilon.internalField();

    const volScalarField& k = db().lookupObject<volScalarField>(kName_);
    const scalarField& kIn = k.internalField();

    const scalarField& nuw =
        lookupPatchField<volScalarField, scalar>(nuName_);

    const scalarField& nutw =
        lookupPatchField<volScalarField, scalar>(nutName_);

    const fvPatchVectorField& Uw =
        lookupPatchField<volVectorField, vector>(UName_);

    const scalarField magGradUw = mag(Uw.snGrad());

    // Get face cells
    const unallocLabelList& fc = patch().faceCells();

    // Averaged G and epsilon
    forAll (nutw, faceI)
    {
        const label faceCellI = fc[faceI];

        // Viscous sublayer thickness
        const scalar yVis = nuw[faceI]*11.225/(Cmu25*sqrt(k[faceCellI]));

        // Height of the wall adjacent volume
        const scalar yn = 2.0*y[faceI];

        if (y[faceI] > yVis)
        {
            GIn[faceCellI] =
                sqr
                (
                    (nutw[faceI] + nuw[faceI])*magGradUw[faceI]
                )
               *log(yn/yVis)
                /(Cmu25*sqrt(kIn[faceCellI])*kappa_*yn);

            epsilonIn[faceCellI] =
                (
                    2.0*nuw[faceI]*kIn[faceCellI]/yVis+pow(kIn[faceCellI], 1.5)
                   *Cmu75*log(yn/yVis)/kappa_
                )/yn;
        }
        else
        {
            GIn[faceCellI] = 0.0;
            epsilonIn[faceCellI] = 2.0*nuw[faceI]*kIn[faceCellI]/sqr(yVis);
        }
    }

    // TODO: perform averaging for cells sharing more than one boundary face

    zeroGradientFvPatchScalarField::updateCoeffs();
}

void kNonEqWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    zeroGradientFvPatchScalarField::write(os);
    this->writeEntry("value", os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "k", "k", kName_);
    writeEntryIfDifferent<word>(os, "epsilon", "epsilon", epsilonName_);
    writeEntryIfDifferent<word>(os, "G", "RASModel::G", GName_);
    writeEntryIfDifferent<word>(os, "nu", "nu", nuName_);
    writeEntryIfDifferent<word>(os, "nut", "nut", nutName_);
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    kNonEqWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
