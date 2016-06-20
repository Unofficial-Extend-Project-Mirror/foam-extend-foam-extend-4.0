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

#include "tractionDisplacementIncrementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "constitutiveModel.H"
#include "solidSolver.H"
#include "unsIncrTotalLagrangianSolid.H"
#include "nonLinearGeometry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tractionDisplacementIncrementFvPatchVectorField::
tractionDisplacementIncrementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    secondOrder_(false)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


tractionDisplacementIncrementFvPatchVectorField::
tractionDisplacementIncrementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_("traction", dict, p.size()),
    pressure_("pressure", dict, p.size()),
    secondOrder_(false)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;

    if (dict.found("secondOrder"))
    {
        secondOrder_ = Switch(dict.lookup("secondOrder"));
        Info << "Second order correction: " << secondOrder_ << endl;
    }

    Info << "Creating traction displacement incr boundary conditions" << endl;
}


tractionDisplacementIncrementFvPatchVectorField::
tractionDisplacementIncrementFvPatchVectorField
(
    const tractionDisplacementIncrementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    traction_(tdpvf.traction_, mapper),
    pressure_(tdpvf.pressure_, mapper),
    secondOrder_(tdpvf.secondOrder_)
{}


tractionDisplacementIncrementFvPatchVectorField::
tractionDisplacementIncrementFvPatchVectorField
(
    const tractionDisplacementIncrementFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_),
    secondOrder_(tdpvf.secondOrder_)
{}


tractionDisplacementIncrementFvPatchVectorField::
tractionDisplacementIncrementFvPatchVectorField
(
    const tractionDisplacementIncrementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_),
    secondOrder_(tdpvf.secondOrder_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionDisplacementIncrementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void tractionDisplacementIncrementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const tractionDisplacementIncrementFvPatchVectorField& dmptf =
        refCast<const tractionDisplacementIncrementFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void tractionDisplacementIncrementFvPatchVectorField::updateCoeffs()
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

    if (isA<solidSolvers::unsIncrTotalLagrangianSolid>(stress))
    {
        const fvsPatchField<tensor>& gradD =
            patch().lookupPatchField<surfaceTensorField, tensor>
            (
                "gradDf"
            );

        nonLinearGeometry::nonLinearType nonLinear =
            nonLinearGeometry::nonLinearNames_.read
            (
                stress.solidProperties().lookup("nonLinear")
            );

        // Switch nonLinear
        // (
        //     stress.solidProperties().lookup("nonLinear")
        // );

        const fvsPatchField<symmTensor>& sigma =
            patch().lookupPatchField<surfaceSymmTensorField, symmTensor>
            (
                "sigmaf"
            );

        vectorField t = traction_;

        // if (nonLinear)
        if (nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN)
        {
            tensorField F = I + gradD.T() + gradDD.T();

            scalarField J = det(F);

            tensorField invF = hinv(F);

            scalarField SoS0 = mag(J*(n & invF));

            vectorField nCurrent = (n & invF);
            nCurrent /= mag(nCurrent);

            // Cauchy traction
            t -= pressure_*nCurrent;

            // 2nd Piola-Kirchhoff traction
            t = (invF & t)*SoS0;
        }
        else
        {
            t -= pressure_*n;
        }

        // 2nd Piola-Kirchhoff traction increment
        vectorField DTraction = t - (n&sigma);

        gradient() =
            DTraction
          - (n & (mu*gradDD.T() - (mu + lambda)*gradDD))
          - n*lambda*tr(gradDD);

        // if (nonLinear)
        if (nonLinear == nonLinearGeometry::TOTAL_LAGRANGIAN)
        {
            gradient() -=
                (n & (mu*(gradDD & gradDD.T())))
              + (n & (mu*(gradDD & gradD.T())))
              + (n & (mu*(gradD & gradDD.T())))
              + 0.5*n*lambda*tr(gradDD & gradDD.T())
              + 0.5*n*lambda*tr(gradDD & gradD.T())
              + 0.5*n*lambda*tr(gradD & gradDD.T());
        }

        if (stress.rheology().plasticityActive())
        {
            gradient() +=
                2*mu
               *(
                    n
                  & stress.rheology()
                   .DEpsilonP().boundaryField()[patch().index()]
                );
        }

        gradient() /= (2.0*mu + lambda);
    }
//     else if  (isA<solidSolvers::unsITLLSSolid>(stress))
//     {
//         const fvsPatchField<tensor>& gradD =
//             patch().lookupPatchField<surfaceTensorField, tensor>
//             (
//                 "gradDf"
//             );

//         const fvsPatchField<tensor>& F =
//             patch().lookupPatchField<surfaceTensorField, tensor>
//             (
//                 "Ff"
//             );

//         const fvsPatchField<vector>& S =
//             patch().lookupPatchField<surfaceVectorField, vector>
//             (
//                 "Sf"
//             );

//         // Full Cauchy traction

//         vectorField t = traction_;

// //         tensorField relF = I + gradDD.T();
// //         tensorField invRelF = hinv(relF);
// //         scalarField relJ = det(relF);

// //         vectorField nCurrent = relJ*(n & invRelF);
// //         nCurrent /= mag(nCurrent);


//          tensorField oldF = I + gradD.T();
// //          scalarField oldJ = det(oldF);
//          tensorField invOldF = hinv(oldF);

// //         tensorField F = I + gradD.T() + gradDD.T();

//         scalarField J = det(F);

// //         tensorField invF = hinv(F);

//         tensorField relF = (F & invOldF);
// //         scalarField relJ = det(relF);
//         tensorField invRelF = hinv(relF);

//         vectorField nCurrent = (invRelF.T() & S);
// //         vectorField nCurrent = (S & (oldF & invF))*(J/oldJ);
//         nCurrent /= mag(nCurrent);

// //         vectorField nCurrent = (n & invF);
// //         nCurrent /= mag(nCurrent);

//         t -= pressure_*nCurrent;


// //         // Kirchhoff traction

// //         vectorField tK = J*t;

// //         // Displacement increment

// //         scalarField K = (lambda + (2.0/3.0)*mu);
// //         symmTensorField b = symm(F & F.T());

// //         tensorField newGradD = gradD + gradDD;

// //         gradient() =
// //             pow(J, 2.0/3.0)*tK
// //           + (1.0/3.0)*mu*tr(b)*n
// //           - 0.5*K*(sqr(J) - 1.0)*pow(J, 2.0/3.0)*n
// //           - mu*n
// //           - (n & (mu*newGradD.T()))
// //           - (n & (mu*(newGradD.T() & newGradD)))
// //           + (n & ((mu + lambda)*newGradD));

// //         gradient() /= (2*mu + lambda);

// //         gradient() -= (n & gradD);

//         // Displacement increment gradient

//         symmTensorField& tau =
//             const_cast<fvsPatchField<symmTensor>&>
//             (
//                 patch().lookupPatchField<surfaceSymmTensorField, symmTensor>
//                 (
//                     "tauf"
//                 )
//             );

// //         scalarField K = (lambda + (2.0/3.0)*mu);
// //         symmTensorField devB = dev(symm(F & F.T()))*pow(J, -2.0/3.0);
// //         symmTensorField tau = mu*devB + 0.5*K*(sqr(J) - 1.0)*I;

//         vectorField uN = S/mag(S);

//         tensorField uGradDD = (F & invOldF) - I;
//         uGradDD = uGradDD.T();

//         gradient() =
//             t - (nCurrent & tau/J)
// //           + (2.0*mu + lambda)*(n & newGradD);
//           + (2.0*mu + lambda)*(uN & uGradDD);

//         gradient() /= (2*mu + lambda);


//         // Transform gradient to initial configuration

//         uGradDD -= uN*(uN & uGradDD);
//         uGradDD += uN*gradient();

//         tensorField newGradDD = ((I + uGradDD.T()) & oldF) - oldF;
//         newGradDD = newGradDD.T();

//         gradient() = (n & newGradDD);


// //         tau -= symm(transform(nCurrent*nCurrent, tau));
// //         tau += J*symm(nCurrent*t);
//     }
    else
    {
        Info << "tractionDisplacementIncrement boundary condition " <<
            "is not defined for solid solver: " << stress.type() << endl;
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}

void tractionDisplacementIncrementFvPatchVectorField
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

    if (secondOrder_)
    {
        vectorField dDDP = (k&gradDD.patchInternalField());

        vectorField nGradDDP = (n&gradDD.patchInternalField());

        Field<vector>::operator=
        (
            this->patchInternalField() + dDDP
          + 0.5*(gradient() + nGradDDP)
           /this->patch().deltaCoeffs()
        );
    }
    else
    {
//         vectorField nGradDDP = (n&gradDD.patchInternalField());

        Field<vector>::operator=
        (
            this->patchInternalField()
          + (k&gradDD.patchInternalField())
//           + nGradDDP/this->patch().deltaCoeffs()
          + gradient()/this->patch().deltaCoeffs()
        );
    }

//     // Looking up solid solver
//     const solidSolver& stress =
//         this->db().objectRegistry::lookupObject<solidSolver>
//         (
//             "solidProperties"
//         );

//     if (isA<solidSolvers::pRveUnsULLSSolid>(stress))
//     {
//         // Scala traction boundary displacement

//         vectorField Sf = this->patch().Sf();
//         vectorField nf = this->patch().nf();

//         vectorField& DDf = *this;

//         scalarField phi = (DDf.component(1) * Sf.component(1));

//         scalarField weights = mag(phi);

//         if(mag(gSum(weights)) > VSMALL)
//         {
//             weights /= gSum(weights);
//         }

//         phi -= weights*gSum(phi);

//         DDf -= nf*(nf & DDf);
//         DDf += nf*phi/mag(Sf);

//         scalar sum = gSum(DDf & Sf );

//         Info << "sum: " << sum << endl;
//     }

    fvPatchField<vector>::evaluate();
}


// Write
void tractionDisplacementIncrementFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);
    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);
    writeEntry("value", os);

    os.writeKeyword("secondOrder")
        << secondOrder_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    tractionDisplacementIncrementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
