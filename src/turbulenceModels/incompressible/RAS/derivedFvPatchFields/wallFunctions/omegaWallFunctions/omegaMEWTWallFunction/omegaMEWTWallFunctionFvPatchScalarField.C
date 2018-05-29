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

#include "omegaMEWTWallFunctionFvPatchScalarField.H"
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

void omegaMEWTWallFunctionFvPatchScalarField::checkType()
{
    if (!patch().isWall())
    {
        FatalErrorIn("omegaMEWTWallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

omegaMEWTWallFunctionFvPatchScalarField::omegaMEWTWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF),
    pName_("p"),
    UName_("U"),
    kName_("k"),
    GName_("RASModel::G"),
    nuName_("nu"),
    nutName_("nut"),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    beta1_(0.075)
{
    checkType();
}


omegaMEWTWallFunctionFvPatchScalarField::omegaMEWTWallFunctionFvPatchScalarField
(
    const omegaMEWTWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedInternalValueFvPatchField<scalar>(ptf, p, iF, mapper),
    pName_(ptf.pName_),
    UName_(ptf.UName_),
    kName_(ptf.kName_),
    GName_(ptf.GName_),
    nuName_(ptf.nuName_),
    nutName_(ptf.nutName_),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    beta1_(ptf.beta1_)
{
    checkType();
}


omegaMEWTWallFunctionFvPatchScalarField::omegaMEWTWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF, dict),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    kName_(dict.lookupOrDefault<word>("k", "k")),
    GName_(dict.lookupOrDefault<word>("G", "RASModel::G")),
    nuName_(dict.lookupOrDefault<word>("nu", "nu")),
    nutName_(dict.lookupOrDefault<word>("nut", "nut")),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    beta1_(dict.lookupOrDefault<scalar>("beta1", 0.075))
{
    checkType();
}


omegaMEWTWallFunctionFvPatchScalarField::omegaMEWTWallFunctionFvPatchScalarField
(
    const omegaMEWTWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedInternalValueFvPatchField<scalar>(owfpsf),
    pName_(owfpsf.pName_),
    UName_(owfpsf.UName_),
    kName_(owfpsf.kName_),
    GName_(owfpsf.GName_),
    nuName_(owfpsf.nuName_),
    nutName_(owfpsf.nutName_),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    beta1_(owfpsf.beta1_)
{
    checkType();
}


omegaMEWTWallFunctionFvPatchScalarField::omegaMEWTWallFunctionFvPatchScalarField
(
    const omegaMEWTWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(owfpsf, iF),
    pName_(owfpsf.pName_),
    UName_(owfpsf.UName_),
    kName_(owfpsf.kName_),
    GName_(owfpsf.GName_),
    nuName_(owfpsf.nuName_),
    nutName_(owfpsf.nutName_),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    beta1_(owfpsf.beta1_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void omegaMEWTWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // If G field is not present, execute zero gradient evaluation
    // HJ, 20/Mar/2011
    if (!db().foundObject<volScalarField>(GName_))
    {
        InfoIn("void omegaMEWTWallFunctionFvPatchScalarField::updateCoeffs()")
            << "Cannot access " << GName_ << " field for patch "
            << patch().name() << ".  Evaluating as zeroGradient"
            << endl;

        fvPatchScalarField::updateCoeffs();
        zeroGradientFvPatchScalarField::evaluate();

        return;
    }

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalarField& y = rasModel.y()[patch().index()];

    volScalarField& G = const_cast<volScalarField&>
        (db().lookupObject<volScalarField>(GName_));

    // Note: omega is now a refValue and set in fixedInternalValueFvPatchField
    // HJ, 3/Aug/2011
    scalarField& omega = refValue();

    const scalarField& nuw =
        lookupPatchField<volScalarField, scalar>(nuName_);

    const scalarField& nutw =
        lookupPatchField<volScalarField, scalar>(nutName_);

    // Velocity field at the patch
    const fvPatchVectorField& Uw =
        lookupPatchField<volVectorField, vector>(UName_);

    // Patch normals
    const vectorField n = patch().nf();

    // Velocity patch internal field
    const vectorField UwIn = Uw.patchInternalField();

    // Velocity vector tangential to the wall
    const vectorField UwInTang = UwIn - (UwIn & n)*n;
    const scalarField magUwInTang = mag(UwInTang);

    // Calculate tangential direction for patch cells
    const vectorField tDir = UwInTang/magUwInTang;

    // Magnitude of the surface normal gradient velocity
    const scalarField magGradUw = mag(Uw.snGrad());


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
        Info<< "Field pressureGradient not found. Neglecting pressure gradient "
            << "effects for wall functions at patch: " << patch().name()
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
        Info<< "Field velocityConvection not found. Neglecting convection "
            << "effects for wall functions at patch: " << patch().name()
            << endl;
    }


    // Get face cells
    const unallocLabelList& fc = patch().faceCells();

    // Set omega and G
    forAll(nutw, faceI)
    {
        const label faceCellI = fc[faceI];

        // Tangential stress at the wall
        const scalar tauw = (nuw[faceI] + nutw[faceI])*magGradUw[faceI];

        const scalar yPlus = sqrt(tauw)*y[faceI]/nuw[faceI];

        // Velocity gradient for viscous sublayer
        const scalar dudyVis= magGradUw[faceI];

        // Velocity gradient for log layer
        const scalar dudyLog =
            sqrt
            (
                max
                (
                    SMALL,
                    (gradpTang[faceI] + convectionTang[faceI])*y[faceI] +tauw
                )
            )/(kappa_*y[faceI]);

        // Kader blending for velocity gradient
        const scalar gamma = -0.01*pow(yPlus, 4)/(1.0 + 5.0*yPlus);
        const scalar dudy = dudyVis*exp(gamma) + dudyLog*exp(1.0/gamma);

        const scalar omegaVis = 6.0*nuw[faceI]/(beta1_*sqr(y[faceI]));
        const scalar omegaLog = dudyLog/(sqrt(Cmu_));

        // Menter blend for omega
        omega[faceI] = sqrt(sqr(omegaVis) + sqr(omegaLog));

        // Production term
        G[faceCellI]= tauw*dudy;
    }

    // TODO: perform averaging for cells sharing more than one boundary face

    fixedInternalValueFvPatchField<scalar>::updateCoeffs();
}

void omegaMEWTWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedInternalValueFvPatchField<scalar>::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntryIfDifferent<word>(os, "k", "k", kName_);
    writeEntryIfDifferent<word>(os, "G", "RASModel::G", GName_);
    writeEntryIfDifferent<word>(os, "nu", "nu", nuName_);
    writeEntryIfDifferent<word>(os, "nut", "nut", nutName_);
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("beta1") << beta1_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    omegaMEWTWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
