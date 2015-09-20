/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "immersedBoundaryEpsilonWallFunctionFvPatchScalarField.H"
#include "immersedBoundaryVelocityWallFunctionFvPatchVectorField.H"
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

immersedBoundaryEpsilonWallFunctionFvPatchScalarField::
immersedBoundaryEpsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(p, iF),
    UName_("U"),
    kName_("k"),
    GName_("RASModel::G"),
    nuName_("nu"),
    nutName_("nut"),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8)
{}


immersedBoundaryEpsilonWallFunctionFvPatchScalarField::
immersedBoundaryEpsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    kName_(dict.lookupOrDefault<word>("k", "k")),
    GName_(dict.lookupOrDefault<word>("G", "RASModel::G")),
    nuName_(dict.lookupOrDefault<word>("nu", "nu")),
    nutName_(dict.lookupOrDefault<word>("nut", "nut")),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8))
{}


immersedBoundaryEpsilonWallFunctionFvPatchScalarField::
immersedBoundaryEpsilonWallFunctionFvPatchScalarField
(
    const immersedBoundaryEpsilonWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    kName_(ptf.kName_),
    GName_(ptf.GName_),
    nuName_(ptf.nuName_),
    nutName_(ptf.nutName_),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{}


immersedBoundaryEpsilonWallFunctionFvPatchScalarField::
immersedBoundaryEpsilonWallFunctionFvPatchScalarField
(
    const immersedBoundaryEpsilonWallFunctionFvPatchScalarField& ewfpsf
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(ewfpsf),
    UName_(ewfpsf.UName_),
    kName_(ewfpsf.kName_),
    GName_(ewfpsf.GName_),
    nuName_(ewfpsf.nuName_),
    nutName_(ewfpsf.nutName_),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_),
    E_(ewfpsf.E_)
{}


immersedBoundaryEpsilonWallFunctionFvPatchScalarField::
immersedBoundaryEpsilonWallFunctionFvPatchScalarField
(
    const immersedBoundaryEpsilonWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    immersedBoundaryWallFunctionFvPatchScalarField(ewfpsf, iF),
    UName_(ewfpsf.UName_),
    kName_(ewfpsf.kName_),
    GName_(ewfpsf.GName_),
    nuName_(ewfpsf.nuName_),
    nutName_(ewfpsf.nutName_),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_),
    E_(ewfpsf.E_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void immersedBoundaryEpsilonWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // If G field is not present, execute zero gradient evaluation
    // HJ, 20/Mar/2011
    if (!db().foundObject<volScalarField>(GName_))
    {
        InfoIn
        (
            "void immersedBoundaryEpsilonWallFunctionFvPatchScalarField::"
            "updateCoeffs()"
        )   << "Cannot access " << GName_ << " field.  for patch "
            << patch().name() << ".  Evaluating as regular immersed boundary"
            << endl;

        immersedBoundaryWallFunctionFvPatchScalarField::evaluate();

        return;
    }

    const vectorField& n = ibPatch().ibNormals();

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalar yPlusLam = rasModel.yPlusLam(kappa_, E_);

    const scalar Cmu25 = pow(Cmu_, 0.25);
    const scalar Cmu50 = sqrt(Cmu_);
    const scalar Cmu75 = pow(Cmu_, 0.75);

    volScalarField& G = const_cast<volScalarField&>
        (db().lookupObject<volScalarField>(GName_));

    // Grab values of other fields required for wall functions

    // Velocity
    const fvPatchVectorField& Uwg =
      patch().lookupPatchField<volVectorField, vector>(UName_);
    const immersedBoundaryVelocityWallFunctionFvPatchVectorField& Uw =
        refCast<const immersedBoundaryVelocityWallFunctionFvPatchVectorField>
        (
            Uwg
        );

    // Calculate tangential component, taking into account wall velocity
    const vectorField UtanOld =
        (I - sqr(n)) & (Uw.ibSamplingPointValue() - Uw.ibValue());
    const scalarField magUtanOld = mag(UtanOld);

    // Tangential velocity component
    scalarField& UTangentialNew = Uw.wallTangentialValue();

    // Wall shear stress
    vectorField& tauWall = Uw.tauWall();

    // Turbulence kinetic energy
    const fvPatchScalarField& kg =
        patch().lookupPatchField<volScalarField, scalar>(kName_);
    const immersedBoundaryWallFunctionFvPatchScalarField& kw =
        refCast<const immersedBoundaryWallFunctionFvPatchScalarField>(kg);

    // Current and new values of k at sampling point
    scalarField k = kw.ibSamplingPointValue();
    scalarField& kNew = kw.wallValue();

    // Laminar viscosity
    const fvPatchScalarField& nuwg =
        patch().lookupPatchField<volScalarField, scalar>(nuName_);
    const immersedBoundaryFvPatchScalarField& nuw =
        refCast<const immersedBoundaryFvPatchScalarField>(nuwg);
    scalarField nu = nuw.ibCellValue();

    // Turbulent viscosity
    const fvPatchScalarField& nutwg =
       patch().lookupPatchField<volScalarField, scalar>(nutName_);
    const immersedBoundaryWallFunctionFvPatchScalarField& nutw =
        refCast<const immersedBoundaryWallFunctionFvPatchScalarField>(nutwg);

    // New values of nut
    scalarField nutOld = nutw.ibCellValue();
    scalarField& nutNew = nutw.wallValue();

    const scalarField magGradUw = mag(Uw.ibGrad());

    // Get the IB addressing and distance
    const labelList& ibc = ibPatch().ibCells();

    // Distance to sampling point
    const scalarField& ySample = ibPatch().ibSamplingPointDelta();

    // Distance from wall to IB point
    const scalarField& y = ibPatch().ibDelta();

    // Epsilon: store IB cell values for direct insertion
    scalarField epsilonSample = this->ibSamplingPointValue();

    scalarField& epsilonNew = this->wallValue();

    // Mark values to be fixed
    boolList wf(ibc.size(), false);

    // Calculate yPlus for sample points
    scalarField ypd = Cmu25*ySample*sqrt(k)/nu;

    // Calculate wall function conditions
    forAll (ibc, ibCellI)
    {
        const scalar nuLam = nu[ibCellI];

        // Calculate yPlus from k and laminar viscosity for the IB point
        const scalar yPlusSample = ypd[ibCellI];

        scalar uTau;

        if (yPlusSample > yPlusLam)
        {
            // Calculate uTau from log-law, knowing sampled k and U
            uTau = magUtanOld[ibCellI]*kappa_/log(E_*yPlusSample);
        }
        else
        {
            // Sampling point is in laminar sublayer
            // Bug fix: HJ, 11/Aug/2014
            uTau = yPlusSample;

        }

        // Set wall shear stress
        tauWall[ibCellI] = sqr(uTau)*UtanOld[ibCellI]/(magUtanOld[ibCellI] + SMALL);

        // Calculate yPlus for IB point
//         scalar yPlusIB = uTau*y[ibCellI]/nuLam;
        scalar yPlusIB = yPlusSample*y[ibCellI]/ySample[ibCellI];

        // Calculate wall function data in the immersed boundary point
        if (yPlusIB > yPlusLam)
        {
            // Logarithmic region
            wf[ibCellI] = true;

            scalar nutw = nuLam*(yPlusIB*kappa_/log(E_*yPlusIB) - 1);

            // Fix generation even though it if is not used
            G[ibc[ibCellI]] =
                sqr((nutw + nuLam)*magGradUw[ibCellI])/
                (Cmu25*sqrt(k[ibCellI])*kappa_*y[ibCellI]);

            // Log-Law for tangential velocity
            UTangentialNew[ibCellI] =
                min
                (
                    magUtanOld[ibCellI],
                    uTau/kappa_*log(E_*yPlusIB)
                );

            // Calculate turbulent viscosity
            nutNew[ibCellI] = nutw;

            // Calculate k in the IB cell from G = epsilon
            kNew[ibCellI] = (nutw + nuLam)*magGradUw[ibCellI]/Cmu50;

            // Calculate epsilon from yPlus and set it
            epsilonNew[ibCellI] =
                Cmu75*pow(kNew[ibCellI], 1.5)/(kappa_*y[ibCellI]);
        }
        else
        {
            // Laminar sub-layer
            wf[ibCellI] = false;

            // G is zero
            G[ibc[ibCellI]] = 0;

            // Laminar sub-layer for tangential velocity: uPlus = yPlus
            UTangentialNew[ibCellI] = min(magUtanOld[ibCellI], uTau*yPlusIB);

            // Turbulent viscosity is zero
            nutNew[ibCellI] = SMALL;

            // k is zero gradient: use the sampled value
            kNew[ibCellI] = k[ibCellI];

            // Calculate epsilon from yPlus and set it.
            // Note:  calculating equilibrium epsilon in the sub-layer creates
            //        an unrealistic oscillation: use sampled value
            // HJ, 27/Jul/2012
            epsilonNew[ibCellI] = epsilonSample[ibCellI];
        }
    }

//     Info<< "UTangentialNew " << min(UTangentialNew) << " " << max(UTangentialNew) << endl;
//     Info<< "nutNew " << min(nutNew) << " " << max(nutNew) << endl;
//     Info<< "kNew " << min(kNew) << " " << max(kNew) << endl;
//     Info<< "epsilonNew " << min(epsilonNew) << " " << max(epsilonNew) << endl;

    // Set the fields to calculated wall function values
    Uw.wallMask() = true;
    kw.wallMask() = wf;
    nutw.wallMask() = true;
    this->wallMask() = true;

    // Insert epsilon values into the internal field
    immersedBoundaryWallFunctionFvPatchScalarField::updateCoeffs();
}


void immersedBoundaryEpsilonWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Insert epsilon values into the internal field
    this->setIbCellValues(this->wallValue());

    fvPatchScalarField::evaluate(commsType);
}


void immersedBoundaryEpsilonWallFunctionFvPatchScalarField::
write(Ostream& os) const
{
    immersedBoundaryWallFunctionFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "k", "k", kName_);
    writeEntryIfDifferent<word>(os, "G", "RASModel::G", GName_);
    writeEntryIfDifferent<word>(os, "nu", "nu", nuName_);
    writeEntryIfDifferent<word>(os, "nut", "nut", nutName_);
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    immersedBoundaryEpsilonWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
