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

#include "epsilonNonEqWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
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

void epsilonNonEqWallFunctionFvPatchScalarField::checkType()
{
    if (!this->patch().isWall())
    {
        FatalErrorIn("epsilonNonEqWallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void epsilonNonEqWallFunctionFvPatchScalarField::
writeLocalEntries(Ostream& os) const
{
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "k", "k", kName_);
    writeEntryIfDifferent<word>(os, "G", "RASModel::G", GName_);
    writeEntryIfDifferent<word>(os, "nu", "nu", nuName_);
    writeEntryIfDifferent<word>(os, "nut", "nut", nutName_);
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

epsilonNonEqWallFunctionFvPatchScalarField::
epsilonNonEqWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchScalarField(p, iF),
    UName_("U"),
    kName_("k"),
    nuName_("nu"),
    nutName_("nut"),
    Cmu_(0.09),
    kappa_(0.41)
{
    checkType();
}


epsilonNonEqWallFunctionFvPatchScalarField::
epsilonNonEqWallFunctionFvPatchScalarField
(
    const epsilonNonEqWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedInternalValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    kName_(ptf.kName_),
    nuName_(ptf.nuName_),
    nutName_(ptf.nutName_),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_)
{
    checkType();
}


epsilonNonEqWallFunctionFvPatchScalarField::
epsilonNonEqWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedInternalValueFvPatchScalarField(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    kName_(dict.lookupOrDefault<word>("k", "k")),
    nuName_(dict.lookupOrDefault<word>("nu", "nu")),
    nutName_(dict.lookupOrDefault<word>("nut", "nut")),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41))
{
    checkType();
}


epsilonNonEqWallFunctionFvPatchScalarField::
epsilonNonEqWallFunctionFvPatchScalarField
(
    const epsilonNonEqWallFunctionFvPatchScalarField& ewfpsf
)
:
    fixedInternalValueFvPatchScalarField(ewfpsf),
    UName_(ewfpsf.UName_),
    kName_(ewfpsf.kName_),
    nuName_(ewfpsf.nuName_),
    nutName_(ewfpsf.nutName_),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_)
{
    checkType();
}


epsilonNonEqWallFunctionFvPatchScalarField::
epsilonNonEqWallFunctionFvPatchScalarField
(
    const epsilonNonEqWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchScalarField(ewfpsf, iF),
    UName_(ewfpsf.UName_),
    kName_(ewfpsf.kName_),
    nuName_(ewfpsf.nuName_),
    nutName_(ewfpsf.nutName_),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void epsilonNonEqWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalarField& y = rasModel.y()[patch().index()];

    const scalar Cmu25 = pow(Cmu_, 0.25);
    const scalar Cmu75 = pow(Cmu_, 0.75);

    // Epsilon is now a refValue (derived from fixedInternalValueFvPatchField)
    scalarField& epsilon = refValue();

    const volScalarField& k = db().lookupObject<volScalarField>(kName_);

    const scalarField& nuw =
        lookupPatchField<volScalarField, scalar>(nuName_);

     const scalarField& nutw =
        lookupPatchField<volScalarField, scalar>(nutName_);

    // Get face cells
    const unallocLabelList& fc = patch().faceCells();

    // Set epsilon
    forAll(nutw, faceI)
    {
        const label faceCellI = fc[faceI];

        // Viscous sublayer thickness
        const scalar yVis = nuw[faceI]*11.225/(Cmu25*sqrt(k[faceCellI]));

        if (y[faceI] > yVis)
        {
           epsilon[faceI] = Cmu75*pow(k[faceCellI], 1.5)/(kappa_*y[faceI]);
        }
        else
        {
           epsilon[faceI] = 2.0*nuw[faceI]*k[faceCellI]/sqr(y[faceI]);
        }
    }

    // TODO: perform averaging for cells sharing more than one boundary face

    fixedInternalValueFvPatchScalarField::updateCoeffs();
}


void epsilonNonEqWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedInternalValueFvPatchScalarField::evaluate(commsType);
}


void epsilonNonEqWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedInternalValueFvPatchScalarField::write(os);
    writeLocalEntries(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    epsilonNonEqWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
