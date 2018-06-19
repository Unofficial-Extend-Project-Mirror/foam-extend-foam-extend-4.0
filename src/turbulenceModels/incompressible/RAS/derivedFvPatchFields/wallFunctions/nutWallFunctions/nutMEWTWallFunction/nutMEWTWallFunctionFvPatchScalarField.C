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

#include "nutMEWTWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::debug::tolerancesSwitch
Foam::incompressible::RASModels::
nutMEWTWallFunctionFvPatchScalarField::dimlessAFactorTol_
(
    "dimlessAFactorMEWTTolerance",
    1e-6
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> nutMEWTWallFunctionFvPatchScalarField::calcNut() const
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

    // Needed to calculate yPlus
    const scalarField& eddyVis =
        lookupPatchField<volScalarField, scalar>(nutName_);

    tmp<scalarField> tnutw(new scalarField(patch().size(), SMALL));
    scalarField& nutw = tnutw();

    // Get face cells
    const unallocLabelList& fc = patch().faceCells();

    forAll(nutw, faceI)
    {
        const label faceCellI = fc[faceI];
        const scalar uStar = Cmu25*sqrt(k[faceCellI]);

        // Calculate yPlus
        const scalar yPlus =
            sqrt
            (
                (eddyVis[faceI]+nuw[faceI])*magGradUw[faceI]
            )*
            y[faceI]/
           (nuw[faceI] + SMALL);

        // Relative tangential velocity
        const scalar magUrel = magUwInTang[faceI] - magUwTang[faceI];

        // Dimless A factor
        scalar A = nuw[faceI]*
            (gradpTang[faceI] + convectionTang[faceI])/(pow(uStar, 3) + SMALL);

        // Numerical stabilisation of the A factor
        if (A < SMALL)
        {
            A += dimlessAFactorTol_();
        }

        // Helper variables
        const scalar S1 = sqrt(max(SMALL, 1.0 + A*yPlus));
        const scalar p1 = sqrt(max(SMALL, 1.0 +6.0*A));
        const scalar uPlusT =
            (1.0/kappa_)*log(6.0*E_)
           -(1.0/kappa_)*(2.0*p1 + log(mag(p1 - 1.0)) + log(p1 + 1.0));

        // Dimless velocity in log layer
        const scalar uLogPlus =
            (1.0/kappa_)*(2.0*S1 + log(mag(S1 - 1.0)) + log(S1 + 1.0)) + uPlusT;

        // Friction velocity in viscous sublayer
        const scalar uTauVis = sqrt(nuw[faceI]*magUrel/(y[faceI] + SMALL));

        // Friction velocity in log layer
        const scalar uTauLog = magUrel/(uLogPlus + SMALL);

        // Kader blending for friction velocity
        const scalar gamma = -0.01*pow(yPlus, 4)/(1.0 + 5.0*yPlus);
        const scalar uTau = uTauVis*exp(gamma) + uTauLog*exp(1.0/gamma);

        // Need to limit nutw for stability reasons since at some point, yPlus
        // is 0 and Kader blending is not defined
        nutw[faceI] =
            max
            (
                SMALL,
                sqr(uTau)/(magGradUw[faceI] + SMALL) - nuw[faceI]
            );
    }

    return tnutw;
}


void nutMEWTWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    nutkWallFunctionFvPatchScalarField::writeLocalEntries(os);

    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntryIfDifferent<word>(os, "nut", "nut", nutName_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutMEWTWallFunctionFvPatchScalarField::nutMEWTWallFunctionFvPatchScalarField
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


nutMEWTWallFunctionFvPatchScalarField::nutMEWTWallFunctionFvPatchScalarField
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


nutMEWTWallFunctionFvPatchScalarField::nutMEWTWallFunctionFvPatchScalarField
(
    const nutMEWTWallFunctionFvPatchScalarField& ptf,
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


nutMEWTWallFunctionFvPatchScalarField::nutMEWTWallFunctionFvPatchScalarField
(
    const nutMEWTWallFunctionFvPatchScalarField& wfpsf
)
:
    nutkWallFunctionFvPatchScalarField(wfpsf),
    UName_(wfpsf.UName_),
    pName_(wfpsf.pName_),
    nutName_(wfpsf.nutName_)
{}


nutMEWTWallFunctionFvPatchScalarField::nutMEWTWallFunctionFvPatchScalarField
(
    const nutMEWTWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(wfpsf, iF),
    UName_(wfpsf.UName_),
    pName_(wfpsf.pName_),
    nutName_(wfpsf.nutName_)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, nutMEWTWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
