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

#include "nutCWTWallFunctionFvPatchScalarField.H"
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void nutCWTWallFunctionFvPatchScalarField::checkType()
{
    if (!patch().isWall())
    {
        FatalErrorIn("nutCWTWallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


scalar nutCWTWallFunctionFvPatchScalarField::calcYPlusLam
(
    const scalar kappa,
    const scalar E
) const
{
    scalar ypl = 11.0;

    for (int i = 0; i < 10; i++)
    {
        ypl = log(E*ypl)/kappa;
    }

    return ypl;
}


tmp<scalarField> nutCWTWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchI = patch().index();

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalarField& y = rasModel.y()[patchI];
    const tmp<volScalarField> tk = rasModel.k();
    const volScalarField& k = tk();
    const scalarField& nuw = rasModel.nu().boundaryField()[patchI];

    const scalar Cmu25 = pow(Cmu_, 0.25);

    // Get normals
    const vectorField n = patch().nf();

    // Patch velocity field at this wall
    const fvPatchVectorField& Uw =
          lookupPatchField<volVectorField, vector>(UName_);

    const scalarField magGradUw = mag(Uw.snGrad());
    const vectorField UwIn = Uw.patchInternalField();

    // Patch internal velocity field tangential to the wall
    const vectorField UwInTang = UwIn - (UwIn & n)*n;
    const scalarField magUwInTang = mag(UwInTang);

    // Calculate tangential direction for patch cells
    const vectorField tDir = UwInTang/magUwInTang;

    // Wall-velocity vector field tangential to the wall
    const vectorField UwTang = Uw - (Uw & n)*n;
    const scalarField magUwTang = mag(UwTang);


    // Pressure terms
    const volScalarField& p =
        this->dimensionedInternalField().mesh().lookupObject
        <
            volScalarField
        >(pName_);

    // Pressure gradient
    const volVectorField gradp = fvc::grad(p);

    // Pressure gradient in wall adjacent cell
    const vectorField gradPIn =
        gradp.boundaryField()[this->patch().index()].patchInternalField();

    // Pressure gradient projected on the wall parallel velocity
    const scalarField gradpTang= gradPIn & tDir;


    // Convective terms
    const volVectorField& U =
        this->dimensionedInternalField().mesh().lookupObject
        <
            volVectorField
        >(UName_);

    const surfaceScalarField& phi =
        this->dimensionedInternalField().mesh().lookupObject
        <
            surfaceScalarField
        >("phi");

    const volVectorField convection = fvc::div(phi, U);

    const vectorField convectionIn =
        convection.boundaryField()[this->patch().index()].patchInternalField();

    // Convection term projected on the wall parallel velocity
    const scalarField convectionTang = convectionIn & tDir;

    tmp<scalarField> tnutw(new scalarField(patch().size()));
    scalarField& nutw = tnutw();

    // Get face cells
    const unallocLabelList& fc = patch().faceCells();

    forAll(nutw, faceI)
    {
        const label faceCellI = fc[faceI];
        const scalar uStar = Cmu25*sqrt(k[faceCellI]);

	    // Note: here yPlus is actually yStar
        const scalar yPlus = uStar*y[faceI]/nuw[faceI];

        // Relative tangential velocity
        const scalar magUrel = magUwInTang[faceI] - magUwTang[faceI];

        const scalar Cu =  convectionTang[faceI] + gradpTang[faceI];
        const scalar Psi = 1.0 - Cu/(kappa_*uStar*magGradUw[faceI]);

        const scalar tauwVis = nuw[faceI]*magGradUw[faceI];
        const scalar tauwLog = kappa_*uStar*magUrel*Psi/log(E_*yPlus);

        // Kader blending
        const scalar gamma = -0.01*pow(yPlus, 4)/(1.0 + 5.0*yPlus);
        const scalar tauw = tauwVis*exp(gamma) + tauwLog*exp(1.0/gamma);

        nutw[faceI] = tauw/magGradUw[faceI] - nuw[faceI];
    }

    return tnutw;
}


void nutCWTWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntryIfDifferent<word>(os, "nut", "nut", nutName_);
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutCWTWallFunctionFvPatchScalarField::nutCWTWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    pName_("p"),
    nutName_("nut"),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    yPlusLam_(calcYPlusLam(kappa_, E_))
{
    checkType();
}


nutCWTWallFunctionFvPatchScalarField::nutCWTWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    nutName_(dict.lookupOrDefault<word>("nut", "nut")),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    yPlusLam_(calcYPlusLam(kappa_, E_))
{
    checkType();
}


nutCWTWallFunctionFvPatchScalarField::nutCWTWallFunctionFvPatchScalarField
(
    const nutCWTWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    pName_(ptf.pName_),
    nutName_(ptf.nutName_),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    yPlusLam_(ptf.yPlusLam_)
{
    checkType();
}


nutCWTWallFunctionFvPatchScalarField::nutCWTWallFunctionFvPatchScalarField
(
    const nutCWTWallFunctionFvPatchScalarField& wfpsf
)
:
    fixedValueFvPatchScalarField(wfpsf),
    UName_(wfpsf.UName_),
    pName_(wfpsf.pName_),
    nutName_(wfpsf.nutName_),
    Cmu_(wfpsf.Cmu_),
    kappa_(wfpsf.kappa_),
    E_(wfpsf.E_),
    yPlusLam_(wfpsf.yPlusLam_)
{
    checkType();
}


nutCWTWallFunctionFvPatchScalarField::nutCWTWallFunctionFvPatchScalarField
(
    const nutCWTWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wfpsf, iF),
    UName_(wfpsf.UName_),
    pName_(wfpsf.pName_),
    nutName_(wfpsf.nutName_),
    Cmu_(wfpsf.Cmu_),
    kappa_(wfpsf.kappa_),
    E_(wfpsf.E_),
    yPlusLam_(wfpsf.yPlusLam_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void nutCWTWallFunctionFvPatchScalarField::updateCoeffs()
{
    operator==(calcNut());

    fixedValueFvPatchScalarField::updateCoeffs();
}


tmp<scalarField> nutCWTWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchI = patch().index();

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalarField& y = rasModel.y()[patchI];

    const tmp<volScalarField> tk = rasModel.k();
    const volScalarField& k = tk();
    const scalarField kwc = k.boundaryField()[patchI].patchInternalField();
    const scalarField& nuw = rasModel.nu().boundaryField()[patchI];

    return pow(Cmu_, 0.25)*y*sqrt(kwc)/nuw;
}


void nutCWTWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, nutCWTWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
