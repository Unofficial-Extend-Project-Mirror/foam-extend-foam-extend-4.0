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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "pRveTractionDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "constitutiveModel.H"
#include "solidSolver.H"
#include "pRveUnsTotalLagrangianSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pRveTractionDisplacementFvPatchVectorField::
pRveTractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


pRveTractionDisplacementFvPatchVectorField::
pRveTractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_("traction", dict, p.size()),
    pressure_("pressure", dict, p.size())
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;

    Info << "Creating traction displacement boundary conditions" << endl;
}


pRveTractionDisplacementFvPatchVectorField::
pRveTractionDisplacementFvPatchVectorField
(
    const pRveTractionDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    traction_(tdpvf.traction_, mapper),
    pressure_(tdpvf.pressure_, mapper)
{}


pRveTractionDisplacementFvPatchVectorField::
pRveTractionDisplacementFvPatchVectorField
(
    const pRveTractionDisplacementFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


pRveTractionDisplacementFvPatchVectorField::
pRveTractionDisplacementFvPatchVectorField
(
    const pRveTractionDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pRveTractionDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void pRveTractionDisplacementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const pRveTractionDisplacementFvPatchVectorField& dmptf =
        refCast<const pRveTractionDisplacementFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void pRveTractionDisplacementFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Looking up solid solver
    const solidSolver& stress =
        this->db().objectRegistry::lookupObject<solidSolver>
        (
            "solidProperties"
        );

    Switch nonLinear
    (
        stress.solidProperties().lookup("nonLinear")
    );

    Switch enforceLinear
    (
        stress.solidProperties().lookup("enforceLinear")
    );

    bool thermalStress = stress.thermalStress();

    word DName = this->dimensionedInternalField().name();

    const fvsPatchField<tensor>& gradDf =
        patch().lookupPatchField<surfaceTensorField, tensor>
        (
            "grad" + DName + "f"
        );

    const fvsPatchField<scalar>& mu =
        patch().lookupPatchField<surfaceScalarField, scalar>
        (
            "muf"
        );

    const fvsPatchField<scalar>& lambda =
        patch().lookupPatchField<surfaceScalarField, scalar>
        (
            "lambdaf"
        );

    vectorField n = patch().nf();

    vectorField t = traction_;

    // Get average deformation gradient from the solver
    const solidSolvers::pRveUnsTotalLagrangianSolid& pRveStress =
        refCast<const solidSolvers::pRveUnsTotalLagrangianSolid>(stress);

    const tensor& avgF = pRveStress.avgDeformationGradient();

    // Calc average strain
    symmTensor avgE = 0.5*symm(avgF + avgF.T()) - I;
    if (nonLinear && !enforceLinear)
    {
        avgE = (symm(avgF.T()&avgF) - I);
    }


    {
        if (nonLinear && !enforceLinear)
        {
            tensorField F = I + gradDf + (avgF - I).T();

            scalarField J = det(F);

            tensorField invF = hinv(F);

            scalarField SoS0 = mag(J*(invF & n));

            vectorField nCurrent = (invF & n);
            nCurrent /= mag(nCurrent);

            t -= pressure_*nCurrent;

            t = (t & invF)*SoS0;
        }
        else
        {
            t -= pressure_*n;
        }

        symmTensorField avgSigma = 2*mu*avgE + lambda*tr(avgE)*I;
        vectorField avgTraction = (n&avgSigma);

        t -= avgTraction;

        gradient() =
            t
          - (n & (mu*gradDf.T() - (mu + lambda)*gradDf))
          - n*lambda*tr(gradDf);

        if (nonLinear && !enforceLinear)
        {
            gradient() -=
                (n & (mu*(gradDf & gradDf.T())))
              + 0.5*n*lambda*tr(gradDf & gradDf.T());
        }

        if (thermalStress)
        {
            const fvPatchField<scalar>& DT =
                patch().lookupPatchField<volScalarField, scalar>("DT");

            const fvsPatchField<scalar>& threeK =
                patch().lookupPatchField<surfaceScalarField, scalar>
                (
                    "threeKf"
                );

            const fvsPatchField<scalar>& alpha =
                patch().lookupPatchField<surfaceScalarField, scalar>
                (
                    "alphaf"
                );

            gradient() += n*threeK*alpha*DT;
        }

        gradient() /= (2.0*mu + lambda);
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}

void pRveTractionDisplacementFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    vectorField n = patch().nf();
    vectorField delta = patch().delta();
    vectorField k = delta - n*(n&delta);

    word UName = this->dimensionedInternalField().name();

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + UName + ")"
        );

    Field<vector>::operator=
    (
        this->patchInternalField()
      + (k&gradU.patchInternalField())
      + gradient()/this->patch().deltaCoeffs()
    );

//     vectorField dUP = (k&gradU.patchInternalField());

//     vectorField nGradUP = (n&gradU.patchInternalField());

//     Field<vector>::operator=
//     (
//         this->patchInternalField() + dUP
//       + 0.5*(gradient() + nGradUP)
//        /this->patch().deltaCoeffs()
//     );

    fvPatchField<vector>::evaluate();
}


// Write
void pRveTractionDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);
    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);
    writeEntry("value", os);
}


// void pRveTractionDisplacementFvPatchVectorField::operator=
// (
//     const fvPatchField<vector>& ptf
// )
// {
//     fvPatchField<vector>::operator=(ptf);

//     gradient() = fvPatchField::snGrad();

//     if (ptf.type() == pRveTractionDisplacementFvPatchVectorField::typeName)
//     {
//         const pRveTractionDisplacementFvPatchVectorField& dmptf =
//             refCast<const pRveTractionDisplacementFvPatchVectorField>(ptf);

//         gradient() = dmptf.gradient();
//     }
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    pRveTractionDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
