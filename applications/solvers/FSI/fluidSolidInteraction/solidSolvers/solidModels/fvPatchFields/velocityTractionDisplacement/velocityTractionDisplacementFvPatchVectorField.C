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

#include "velocityTractionDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "constitutiveModel.H"
#include "solidSolver.H"
// #include "unsTLPressureDisplacementSolid.H"
#include "unsTotalLagrangianSolid.H"
// #include "unsTLLSSolid.H"
// #include "unsTLViscoelasticSolid.H"

#include "tractionBoundaryGradient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

velocityTractionDisplacementFvPatchVectorField::
velocityTractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    robinFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    velocity_(p.size(), vector::zero),
    normal_(p.size(), vector::zero),
    secondOrder_(false),
    limitCoeff_(1.0)
{
    fvPatchVectorField::operator=(patchInternalField());
//     gradient() = vector::zero;
}


velocityTractionDisplacementFvPatchVectorField::
velocityTractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    robinFvPatchVectorField(p, iF),
    traction_("traction", dict, p.size()),
    pressure_("pressure", dict, p.size()),
    velocity_("velocity", dict, p.size()),
    normal_(p.size(), vector::zero),
    secondOrder_(false),
    limitCoeff_(1.0)
{
    fvPatchVectorField::operator=(patchInternalField());
//     gradient() = vector::zero;

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

    Info << "Creating velocity traction displacement boundary conditions"
        << endl;
}


velocityTractionDisplacementFvPatchVectorField::
velocityTractionDisplacementFvPatchVectorField
(
    const velocityTractionDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    robinFvPatchVectorField(tdpvf, p, iF, mapper),
    traction_(tdpvf.traction_, mapper),
    pressure_(tdpvf.pressure_, mapper),
    velocity_(tdpvf.velocity_, mapper),
    normal_(tdpvf.normal_, mapper),
    secondOrder_(tdpvf.secondOrder_),
    limitCoeff_(tdpvf.limitCoeff_)
{}


velocityTractionDisplacementFvPatchVectorField::
velocityTractionDisplacementFvPatchVectorField
(
    const velocityTractionDisplacementFvPatchVectorField& tdpvf
)
:
    robinFvPatchVectorField(tdpvf),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_),
    velocity_(tdpvf.velocity_),
    normal_(tdpvf.normal_),
    secondOrder_(tdpvf.secondOrder_),
    limitCoeff_(tdpvf.limitCoeff_)
{}


velocityTractionDisplacementFvPatchVectorField::
velocityTractionDisplacementFvPatchVectorField
(
    const velocityTractionDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    robinFvPatchVectorField(tdpvf, iF),
    traction_(tdpvf.traction_),
    pressure_(tdpvf.pressure_),
    velocity_(tdpvf.velocity_),
    normal_(tdpvf.normal_),
    secondOrder_(tdpvf.secondOrder_),
    limitCoeff_(tdpvf.limitCoeff_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void velocityTractionDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    robinFvPatchVectorField::autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
    velocity_.autoMap(m);
    normal_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void velocityTractionDisplacementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    robinFvPatchVectorField::rmap(ptf, addr);

    const velocityTractionDisplacementFvPatchVectorField& dmptf =
        refCast<const velocityTractionDisplacementFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
    velocity_.rmap(dmptf.velocity_, addr);
    normal_.rmap(dmptf.normal_, addr);
}


// Update the coefficients associated with the patch field
void velocityTractionDisplacementFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

//     Info << "velocityTractionDisplacementFvPatchVectorField::updateCoeffs"
//         << endl;

    word fieldName = dimensionedInternalField().name();

    word totalFieldName("D");
    if (fieldName == "DU")
    {
        totalFieldName = word("U");
    }

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    scalar dt = mesh.time().deltaT().value();

//     // Looking up solid solver
//     const solidSolver& stress =
//         this->db().objectRegistry::lookupObject<solidSolver>
//         (
//             "solidProperties"
//         );

//     scalar rho = stress.rheology().rho()()[0];
//     scalar mu = stress.rheology().mu()()[0];
//     scalar lambda = stress.rheology().lambda()()[0];
//     scalar ap = sqrt((lambda+2*mu)/rho);
//     scalar h = ap*dt;

    const volVectorField& D =
        mesh.lookupObject<volVectorField>(fieldName);

    const fvsPatchField<tensor>& gradDf =
        patch().lookupPatchField<surfaceTensorField, tensor>
        (
            "grad" + fieldName + "f"
        );

    const fvPatchField<scalar>& rho =
        patch().lookupPatchField<volScalarField, scalar>
        (
            "rho"
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

    const vectorField& pD =
        D.boundaryField()[patch().index()];
    // const vectorField& pOldD =
    //     D.oldTime().boundaryField()[patch().index()];
    const vectorField& pOldOldD =
        D.oldTime().oldTime().boundaryField()[patch().index()];

    scalarField ap = sqrt((lambda+2*mu)/rho);
    scalarField h = ap*dt;

    // Ref configuration normal
    vectorField n = this->patch().nf();

    vectorField nq =
        (n & (mu*gradDf.T() - (mu + lambda)*gradDf))
      + n*lambda*tr(gradDf);

//     if (true)
//     {
//         nq += (n & (mu*(gradDf & gradDf.T())))
//           + 0.5*n*lambda*tr(gradDf & gradDf.T());
//     }


//     this->coeff0() = 0;
//     this->coeff1() = (2*mu+lambda);
//     this->rhs() =
//         (
//             traction_
//           + rho*h*transform(sqr(normal_), dt*velocity_ - (pOldD-pOldOldD))
//            /sqr(dt)
//           - nq
//         );


    this->coeff0() = 1.0;
    this->coeff1() = (2*mu+lambda)*sqr(dt)/(rho*h);
    this->rhs() =
        (
            traction_
          + rho*h*transform(sqr(normal_), pOldOldD + 2*dt*velocity_)/sqr(dt)
          - nq
        )
       *sqr(dt)/(rho*h)
      + transform(I-sqr(normal_), pD);




//     gradient() =
//         tractionBoundaryGradient().snGrad
//         (
//             traction_ + tractionCorr,
//             pressure_,
//             fieldName,
//             totalFieldName,
//             patch(),
//             bool((fieldName == "DU") || (fieldName == "DD"))
//         );

//     scalarField dn = 1.0/this->patch().deltaCoeffs();

//     Info << "Dcoeff " << max(this->coeff()+dn)
//         << ", " << average(this->coeff()+dn) << endl;

    robinFvPatchVectorField::updateCoeffs();
}

// void velocityTractionDisplacementFvPatchVectorField::evaluate
// (
//     const Pstream::commsTypes
// )
// {
//     if (!this->updated())
//     {
//         this->updateCoeffs();
//     }

//     vectorField n = patch().nf();
//     vectorField delta = patch().delta();
//     vectorField k = delta - n*(n&delta);

//     word UName = this->dimensionedInternalField().name();

//     const fvPatchField<tensor>& gradU =
//         patch().lookupPatchField<volTensorField, tensor>
//         (
//             "grad(" + UName + ")"
//         );

//     if (secondOrder_)
//     {
//         vectorField dUP = (k&gradU.patchInternalField());
//         vectorField nGradUP = (n&gradU.patchInternalField());

//         Field<vector>::operator=
//         (
//             this->patchInternalField() + dUP
//           + 0.5*(gradient() + nGradUP)
//            /this->patch().deltaCoeffs()
//         );
//     }
//     else
//     {
//         Field<vector>::operator=
//         (
//             this->patchInternalField()
//           + (k&gradU.patchInternalField())
//           + gradient()/this->patch().deltaCoeffs()
//         );
//     }

//     fvPatchField<vector>::evaluate();
// }


// Write
void velocityTractionDisplacementFvPatchVectorField::write(Ostream& os) const
{
    robinFvPatchVectorField::write(os);
    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);
    velocity_.writeEntry("velocity", os);
    normal_.writeEntry("normal", os);

    os.writeKeyword("secondOrder")
        << secondOrder_ << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    velocityTractionDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
