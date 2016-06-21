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

#include "tractionDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "constitutiveModel.H"
#include "solidSolver.H"
#include "unsTotalLagrangianSolid.H"

#include "tractionBoundaryGradient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    secondOrder_(false),
    limitCoeff_(1.0)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_("traction", dict, p.size()),
    pressure_("pressure", dict, p.size()),
    secondOrder_(false),
    limitCoeff_(1.0)
{
    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
        gradient() = vector::zero;
    }

    if (dict.found("secondOrder"))
    {
        secondOrder_ = Switch(dict.lookup("secondOrder"));
        Info << "Second order correction: " << secondOrder_ << endl;
    }

    if (dict.found("limitCoeff"))
    {
        limitCoeff_ =
            scalar(readScalar(dict.lookup("limitCoeff")));
        Info << "Limiter coefficient: " << limitCoeff_ << endl;
    }

    Info << "Creating traction displacement boundary conditions" << endl;
}


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    traction_(tdpvf.traction_, mapper),
    pressure_(tdpvf.pressure_, mapper),
    secondOrder_(tdpvf.secondOrder_),
    limitCoeff_(tdpvf.limitCoeff_)
{}


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_),
    secondOrder_(tdpvf.secondOrder_),
    limitCoeff_(tdpvf.limitCoeff_)
{}


tractionDisplacementFvPatchVectorField::
tractionDisplacementFvPatchVectorField
(
    const tractionDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_),
    secondOrder_(tdpvf.secondOrder_),
    limitCoeff_(tdpvf.limitCoeff_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void tractionDisplacementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const tractionDisplacementFvPatchVectorField& dmptf =
        refCast<const tractionDisplacementFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void tractionDisplacementFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Info << "tractionDisplacementFvPatchVectorField::updateCoeffs()" << endl;

    word fieldName = dimensionedInternalField().name();

    word totalFieldName("D");
    if (fieldName == "DU")
    {
        totalFieldName = word("U");
    }

    gradient() =
        tractionBoundaryGradient().snGrad
        (
            traction_,
            pressure_,
            fieldName,
            totalFieldName,
            patch(),
            bool((fieldName == "DU") || (fieldName == "DD"))
        );



    // // Looking up solid solver
    // const solidSolver& stress =
    //     this->db().objectRegistry::lookupObject<solidSolver>
    //     (
    //         "solidProperties"
    //     );

    // word DName = this->dimensionedInternalField().name();

    // const fvsPatchField<tensor>& gradDf =
    //     patch().lookupPatchField<surfaceTensorField, tensor>
    //     (
    //         "grad" + DName + "f"
    //     );

    // const fvsPatchField<scalar>& mu =
    //     patch().lookupPatchField<surfaceScalarField, scalar>
    //     (
    //         "muf"
    //     );

    // const fvsPatchField<scalar>& lambda =
    //     patch().lookupPatchField<surfaceScalarField, scalar>
    //     (
    //         "lambdaf"
    //     );

    // vectorField n = patch().nf();

    // vectorField t = traction_;

    // if (isA<solidSolvers::unsTotalLagrangianSolid>(stress))
    // {
    //     // gradient() =
    //     //     tractionBoundaryGradient()
    //     //     (
    //     //         traction_,
    //     //         pressure_,
    //     //         word(dimensionedInternalField().name()),
    //     //         patch(),
    //     //         false
    //     //     )();

    //     // Switch nonLinear
    //     // (
    //     //     stress.solidProperties().lookup("nonLinear")
    //     // );

    //     // Switch enforceLinear
    //     // (
    //     //     stress.solidProperties().lookup("enforceLinear")
    //     // );

    //   Switch nonLinear(true);
    //   Switch enforceLinear(false);

    //     if (nonLinear && !enforceLinear)
    //     {
    //         tensorField F = I + gradDf.T();

    //         scalarField J = det(F);

    //         tensorField invF = hinv(F);

    //         scalarField SoS0 = mag(J*(n & invF));

    //         vectorField nCurrent = (n & invF);
    //         nCurrent /= mag(nCurrent);

    //         t -= pressure_*nCurrent;

    //         t = (invF & t)*SoS0;
    //     }
    //     else
    //     {
    //         t -= pressure_*n;
    //     }

    //     gradient() =
    //         t
    //       - (n & (mu*gradDf.T() - (mu + lambda)*gradDf))
    //       - n*lambda*tr(gradDf);

    //     if (nonLinear && !enforceLinear)
    //     {
    //         gradient() -=
    //             (n & (mu*(gradDf & gradDf.T())))
    //           + 0.5*n*lambda*tr(gradDf & gradDf.T());
    //     }

    //     bool thermalStress = stress.thermalStress();

    //     if (thermalStress)
    //     {
    //         const fvsPatchField<scalar>& DT =
    //             patch().lookupPatchField<surfaceScalarField, scalar>("DTf");

    //         const fvsPatchField<scalar>& threeK =
    //             patch().lookupPatchField<surfaceScalarField, scalar>
    //             (
    //                 "threeKf"
    //             );

    //         const fvsPatchField<scalar>& alpha =
    //             patch().lookupPatchField<surfaceScalarField, scalar>
    //             (
    //                 "alphaf"
    //             );

    //         gradient() += n*threeK*alpha*DT;
    //     }

    //     gradient() /= (2.0*mu + lambda);
    // }
//     else if (isA<solidSolvers::unsTLLSSolid>(stress))
//     {
//         // Full Cauchy traction

//         vectorField t = traction_;

//         tensorField F = I + gradDf.T();

//         scalarField J = det(F);

//         tensorField invF = hinv(F);

//         vectorField nCurrent = (n & invF);
//         nCurrent /= mag(nCurrent);

//         t -= pressure_*nCurrent;


//         // Displacement gradient

// //         vectorField oldGradient = gradient();

//         scalarField K = (lambda + (2.0/3.0)*mu);
//         symmTensorField devB = dev(symm(F & F.T()))*pow(J, -2.0/3.0);
//         symmTensorField tau = mu*devB + 0.5*K*(sqr(J)-1.0)*I;

//         gradient() =
//             t - (nCurrent & tau/J)
//           + (2.0*mu + lambda)*(n & gradDf);

// //         symmTensorField b = symm(F&F.T());
// //         gradient() =
// //             pow(J, 2.0/3.0)*tK
// //           + (1.0/3.0)*mu*tr(b)*n
// //           - 0.5*K*(sqr(J) - 1.0)*pow(J, 2.0/3.0)*n
// //           - mu*n
// //           - (n & (mu*gradDf.T()))
// //           - (n & (mu*(gradDf.T() & gradDf)))
// //           + (n & ((mu + lambda)*gradDf));

//         gradient() /= (2*mu + lambda);

// //         gradient() =
// //             oldGradient
// //           + limitCoeff_*(gradient() - oldGradient);

// //         symmTensorField& tauf =
// //             const_cast<fvsPatchField<symmTensor>&>
// //             (
// //                 patch().lookupPatchField<surfaceSymmTensorField, symmTensor>
// //                 (
// //                     "tauf"
// //                 )
// //             );

// //         tauf -= symm(transform(nCurrent*nCurrent, tauf));
// //         tauf += J*symm(nCurrent*t);
//     }
//     else if (isA<solidSolvers::unsTLPressureDisplacementSolid>(stress))
//     {
//         Switch nonLinear
//         (
//             stress.solidProperties().lookup("nonLinear")
//         );

//         Switch enforceLinear
//         (
//             stress.solidProperties().lookup("enforceLinear")
//         );

//         if (nonLinear && !enforceLinear)
//         {
//             tensorField F = I + gradDf;

//             scalarField J = det(F);

//             tensorField invF = hinv(F);

//             scalarField SoS0 = mag(J*(invF & n));

//             vectorField nCurrent = (invF & n);
//             nCurrent /= mag(nCurrent);

//             t -= pressure_*nCurrent;

//             t = (t & invF)*SoS0;
//         }
//         else
//         {
//             t -= pressure_*n;
//         }

//         const fvPatchField<scalar>& p =
//             patch().lookupPatchField<volScalarField, scalar>("p");

//         gradient() =
//             t
//           - (n & (mu*gradDf.T() - (mu + lambda)*gradDf))
// //           + n*(2.0/3.0)*mu*tr(gradDf)
//           + n*p;

//         if (nonLinear && !enforceLinear)
//         {
//             gradient() -= (n & (mu*(gradDf & gradDf.T())));
// //               - n*(1.0/3.0)*mu*tr(gradDf & gradDf.T());
//         }

//         bool thermalStress = stress.thermalStress();

//         if (thermalStress)
//         {
//             const fvPatchField<scalar>& DT =
//                 patch().lookupPatchField<volScalarField, scalar>("DT");

//             const fvsPatchField<scalar>& threeK =
//                 patch().lookupPatchField<surfaceScalarField, scalar>
//                 (
//                     "threeKf"
//                 );

//             const fvsPatchField<scalar>& alpha =
//                 patch().lookupPatchField<surfaceScalarField, scalar>
//                 (
//                     "alphaf"
//                 );

//             gradient() += n*threeK*alpha*DT;
//         }

//         gradient() /= (2.0*mu + lambda);
//     }
//     else if (isA<solidSolvers::unsTLViscoelasticSolid>(stress))
//     {
//         gradient() =
//             tractionBoundaryGradient()
//             (
//                 traction_,
//                 pressure_,
//                 word(dimensionedInternalField().name()),
//                 patch(),
//                 false
//             )();

//         const fvsPatchField<symmTensor>& sigmaCorrf =
//             patch().lookupPatchField<surfaceSymmTensorField, symmTensor>
//             (
//                 "sigmaCorrf"
//             );

//         Switch nonLinear
//         (
//             stress.solidProperties().lookup("nonLinear")
//         );

//         Switch enforceLinear
//         (
//             stress.solidProperties().lookup("enforceLinear")
//         );

//         if (nonLinear && !enforceLinear)
//         {
//             tensorField F = I + gradDf.T();

//             scalarField J = det(F);

//             tensorField invF = hinv(F);

//             scalarField SoS0 = mag(J*(n & invF));

//             vectorField nCurrent = (n & invF);
//             nCurrent /= mag(nCurrent);

//             t -= pressure_*nCurrent;

//             t = (invF & t)*SoS0;
//         }
//         else
//         {
//             t -= pressure_*n;
//         }

//         gradient() =
//             t
//           - (n & (mu*gradDf.T() - (mu + lambda)*gradDf))
//           - n*lambda*tr(gradDf)
//           - (n & sigmaCorrf);

//         if (nonLinear && !enforceLinear)
//         {
//             gradient() -=
//                 (n & (mu*(gradDf & gradDf.T())))
//               + 0.5*n*lambda*tr(gradDf & gradDf.T());
//         }

//         gradient() /= (2*mu + lambda);
//     }
//     else
//     {
//         Info << "tractionDisplacement boundary condition is not defined "
//             << "for solid solver: " << stress.type() << endl;
//     }

    fixedGradientFvPatchVectorField::updateCoeffs();
}

void tractionDisplacementFvPatchVectorField::evaluate
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

    if (secondOrder_)
    {
        vectorField dUP = (k&gradU.patchInternalField());
        vectorField nGradUP = (n&gradU.patchInternalField());

        Field<vector>::operator=
        (
            this->patchInternalField() + dUP
          + 0.5*(gradient() + nGradUP)
           /this->patch().deltaCoeffs()
        );
    }
    else
    {
//         vectorField uncorrectedSnGrad =
//             (
//                 *this
//               - this->patchInternalField()
//             )
//            *this->patch().deltaCoeffs();

//         vectorField snGradCorrection =
//           - (k&gradU.patchInternalField())
//            *this->patch().deltaCoeffs();

//         scalarField limiter =
//         (
//             min
//             (
//                 limitCoeff_*mag(uncorrectedSnGrad + snGradCorrection)
//                /((1 - limitCoeff_)*mag(snGradCorrection) + SMALL),
//                 1.0
//             )
//         );

//         snGradCorrection *= limiter;

//         Field<vector>::operator=
//         (
//             this->patchInternalField()
//           - snGradCorrection/this->patch().deltaCoeffs()
//           + gradient()/this->patch().deltaCoeffs()
//         );

//         Info << "evaluate traction displacement: "
//             << this->patch().name() << endl;

        Field<vector>::operator=
        (
            this->patchInternalField()
          + (k&gradU.patchInternalField())
          + gradient()/this->patch().deltaCoeffs()
        );
    }

    fvPatchField<vector>::evaluate();
}


// Write
void tractionDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);
    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);
    writeEntry("value", os);

    os.writeKeyword("secondOrder")
        << secondOrder_ << token::END_STATEMENT << nl;
}


// void tractionDisplacementFvPatchVectorField::operator=
// (
//     const fvPatchField<vector>& ptf
// )
// {
//     fvPatchField<vector>::operator=(ptf);

//     gradient() = fvPatchField::snGrad();

//     if (ptf.type() == tractionDisplacementFvPatchVectorField::typeName)
//     {
//         const tractionDisplacementFvPatchVectorField& dmptf =
//             refCast<const tractionDisplacementFvPatchVectorField>(ptf);

//         gradient() = dmptf.gradient();
//     }
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, tractionDisplacementFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
