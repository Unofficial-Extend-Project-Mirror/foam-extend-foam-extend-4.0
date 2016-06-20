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

#include "pRveTractionDisplacementIncrementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "constitutiveModel.H"
#include "solidSolver.H"
#include "pRveUnsIncrTotalLagrangianSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pRveTractionDisplacementIncrementFvPatchVectorField::
pRveTractionDisplacementIncrementFvPatchVectorField
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


pRveTractionDisplacementIncrementFvPatchVectorField::
pRveTractionDisplacementIncrementFvPatchVectorField
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

    Info << "Creating pRve traction displacement incr boundary conditions"
        << endl;
}


pRveTractionDisplacementIncrementFvPatchVectorField::
pRveTractionDisplacementIncrementFvPatchVectorField
(
    const pRveTractionDisplacementIncrementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    traction_(tdpvf.traction_, mapper),
    pressure_(tdpvf.pressure_, mapper)
{}


pRveTractionDisplacementIncrementFvPatchVectorField::
pRveTractionDisplacementIncrementFvPatchVectorField
(
    const pRveTractionDisplacementIncrementFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


pRveTractionDisplacementIncrementFvPatchVectorField::
pRveTractionDisplacementIncrementFvPatchVectorField
(
    const pRveTractionDisplacementIncrementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pRveTractionDisplacementIncrementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void pRveTractionDisplacementIncrementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const pRveTractionDisplacementIncrementFvPatchVectorField& dmptf =
        refCast<const pRveTractionDisplacementIncrementFvPatchVectorField>
        (
            ptf
        );

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void pRveTractionDisplacementIncrementFvPatchVectorField::updateCoeffs()
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

//     Switch nonLinear
//     (
//         stress.solidProperties().lookup("nonLinear")
//     );

//     Switch enforceLinear
//     (
//         stress.solidProperties().lookup("enforceLinear")
//     );

    word DDName = this->dimensionedInternalField().name();

    const fvsPatchField<tensor>& gradDD =
        patch().lookupPatchField<surfaceTensorField, tensor>
        (
            "grad" + DDName + "f"
        );

//     const fvsPatchField<tensor>& gradD =
//         patch().lookupPatchField<surfaceTensorField, tensor>
//         (
//             "gradDf"
//         );

    const fvPatchField<symmTensor>& totSigma =
        patch().lookupPatchField<volSymmTensorField, symmTensor>
        (
            "totSigma"
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

    // Get average deformation gradient from the solver
    const solidSolvers::pRveUnsIncrTotalLagrangianSolid& pRveStress =
        refCast<const solidSolvers::pRveUnsIncrTotalLagrangianSolid>
        (
            stress
        );

    const symmTensor& avgDEpsilon = pRveStress.avgDEpsilon();

    {
        vectorField t = traction_;

//         if (nonLinear && !enforceLinear)
//         {
//             tensorField F = I + gradD + gradDD;

//             scalarField J = det(F);

//             tensorField invF = hinv(F);

//             scalarField SoS0 = mag(J*(invF & n));

//             vectorField nCurrent = (invF & n);
//             nCurrent /= mag(nCurrent);

//             // Cauchy traction
//             t -= pressure_*nCurrent;

//             // 2nd Piola-Kirchhoff traction
//             t = (t & invF)*SoS0;
//         }
//         else
        {
            t -= pressure_*n;
        }

        // Total traction increment
        vectorField DTraction = t - (n&totSigma);

        // Avg traction increment
        symmTensorField avgDSigma =
            2*mu*avgDEpsilon + lambda*tr(avgDEpsilon)*I;
        vectorField avgDTraction = (n&avgDSigma);

//         t -= avgTraction;

        // Preturbation traction increment
        DTraction -= avgDTraction;



        gradient() =
            DTraction
          - (n & (mu*gradDD.T() - (mu + lambda)*gradDD))
          - n*lambda*tr(gradDD);

//         if(nonLinear && !enforceLinear)
//         {
//             gradient() -=
//                 (n & (mu*(gradDD & gradDD.T())))
//               + (n & (mu*(gradDD & gradD.T())))
//               + (n & (mu*(gradD & gradDD.T())))
//               + 0.5*n*lambda*tr(gradDD & gradDD.T())
//               + 0.5*n*lambda*tr(gradDD & gradD.T())
//               + 0.5*n*lambda*tr(gradD & gradDD.T());
//         }

        gradient() /= (2.0*mu + lambda);
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}

void pRveTractionDisplacementIncrementFvPatchVectorField
::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    vectorField n = patch().nf();
    vectorField delta = patch().delta();
    vectorField k = delta - n*(n&delta);

    word DDName = this->dimensionedInternalField().name();

    const fvPatchField<tensor>& gradDD =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + DDName + ")"
        );

    Field<vector>::operator=
    (
        this->patchInternalField()
      + (k&gradDD.patchInternalField())
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
void pRveTractionDisplacementIncrementFvPatchVectorField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchVectorField::write(os);
    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    pRveTractionDisplacementIncrementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
