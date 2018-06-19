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

tmp<scalarField> nutCWTWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchI = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const scalarField& y = turbModel.y()[patchI];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const scalarField& nuw = turbModel.nu().boundaryField()[patchI];

    const scalar Cmu25 = pow025(Cmu_);

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
    const vectorField tDir = UwInTang/(magUwInTang + SMALL);

    // Wall-velocity vector field tangential to the wall
    const vectorField UwTang = Uw - (Uw & n)*n;
    const scalarField magUwTang = mag(UwTang);


    // Pressure effects. Lookup the pressure gradient field from the registry
    // (created with pressureGradient function object)

    // Pressure gradient projected on the wall parallel velocity direction
    scalarField gradpTang(magGradUw.size(), 0.0);

    if
    (
        this->dimensionedInternalField().mesh().foundObject
        <
            volVectorField
        >("pressureGradient")
    )
    {
        const volVectorField& gradP =
            this->dimensionedInternalField().mesh().lookupObject
            <volVectorField>("pressureGradient");

        // Update pressure gradient projected on the wall parallel velocity
        gradpTang =
            gradP.boundaryField()[this->patch().index()].patchInternalField()
          & tDir;
    }
    else
    {
        InfoIn
        (
            "tmp<scalarField>"
            "nutCWTWallFunctionFvPatchScalarField::calcNut() const"
        )   << "Field pressureGradient not found. Neglecting pressure gradient "
            << "effects for wall functions at patch: " << patch().name()
            << nl
            << "If you would like to include pressure gradient effects, set up"
            << " pressureGradient function object."
            << endl;
    }


    // Convective terms. Lookup the convection field from the registry (created
    // with velocityConvection function object)

    // Velocity convection projected on the wall parallel velocity direction
    scalarField convectionTang(magGradUw.size(), 0.0);

    if
    (
        this->dimensionedInternalField().mesh().foundObject
        <
            volVectorField
        >("velocityConvection")
    )
    {
        const volVectorField& convection =
            this->dimensionedInternalField().mesh().lookupObject
            <volVectorField>("velocityConvection");

        // Upate convection term projected on the wall parallel velocity
        convectionTang =
            convection.boundaryField()
            [
                this->patch().index()
            ].patchInternalField()
          & tDir;
    }
    else
    {
        InfoIn
        (
            "tmp<scalarField>"
            "nutCWTWallFunctionFvPatchScalarField::calcNut() const"
        )   << "Field velocityConvection not found. Neglecting convection "
            << "effects for wall functions at patch: " << patch().name()
            << nl
            << "If you would like to include convection effects, set up"
            << " velocityConvection function object."
            << endl;
    }


    tmp<scalarField> tnutw(new scalarField(patch().size(), SMALL));
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

        nutw[faceI] = tauw/(magGradUw[faceI] + SMALL) - nuw[faceI];
    }

    return tnutw;
}


void nutCWTWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    nutkWallFunctionFvPatchScalarField::writeLocalEntries(os);

    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntryIfDifferent<word>(os, "nut", "nut", nutName_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutCWTWallFunctionFvPatchScalarField::nutCWTWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(p, iF),
    UName_("U"),
    pName_("p"),
    nutName_("nut")
{}


nutCWTWallFunctionFvPatchScalarField::nutCWTWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutkWallFunctionFvPatchScalarField(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    nutName_(dict.lookupOrDefault<word>("nut", "nut"))
{}


nutCWTWallFunctionFvPatchScalarField::nutCWTWallFunctionFvPatchScalarField
(
    const nutCWTWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    pName_(ptf.pName_),
    nutName_(ptf.nutName_)
{}


nutCWTWallFunctionFvPatchScalarField::nutCWTWallFunctionFvPatchScalarField
(
    const nutCWTWallFunctionFvPatchScalarField& wfpsf
)
:
    nutkWallFunctionFvPatchScalarField(wfpsf),
    UName_(wfpsf.UName_),
    pName_(wfpsf.pName_),
    nutName_(wfpsf.nutName_)
{}


nutCWTWallFunctionFvPatchScalarField::nutCWTWallFunctionFvPatchScalarField
(
    const nutCWTWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(wfpsf, iF),
    UName_(wfpsf.UName_),
    pName_(wfpsf.pName_),
    nutName_(wfpsf.nutName_)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, nutCWTWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
