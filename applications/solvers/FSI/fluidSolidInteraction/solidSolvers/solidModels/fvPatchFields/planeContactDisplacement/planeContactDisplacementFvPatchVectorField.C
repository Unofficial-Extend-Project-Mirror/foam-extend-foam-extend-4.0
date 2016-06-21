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

Description

\*---------------------------------------------------------------------------*/

#include "planeContactDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
// #include "rheologyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

planeContactDisplacementFvPatchVectorField::
planeContactDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    p0_(vector::zero),
    n0_(vector::zero),
    U0_(0.0)
{
    fvPatchVectorField::operator=(patchInternalField());
    refGrad() = vector::zero;
    refValue() = vector::zero;
}


planeContactDisplacementFvPatchVectorField::
planeContactDisplacementFvPatchVectorField
(
    const planeContactDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(tdpvf, p, iF, mapper),
    p0_(tdpvf.p0_),
    n0_(tdpvf.n0_),
    U0_(tdpvf.U0_)
{}


planeContactDisplacementFvPatchVectorField::
planeContactDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    p0_(dict.lookup("p0")),
    n0_(dict.lookup("n0")),
    U0_(readScalar(dict.lookup("U0")))
{
    this->valueFraction() = symmTensor::zero;
    this->refValue() = vector::zero;
    this->refGrad() = vector::zero;

    // Evaluate
    Field<vector> normalValue =
        transform(this->valueFraction(), this->refValue());

    Field<vector> gradValue =
        this->patchInternalField()
      + this->refGrad()/this->patch().deltaCoeffs();

    Field<vector> transformGradValue =
        transform(I - this->valueFraction(), gradValue);

    Field<vector>::operator=(normalValue + transformGradValue);
}


planeContactDisplacementFvPatchVectorField::
planeContactDisplacementFvPatchVectorField
(
    const planeContactDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(tdpvf, iF),
    p0_(tdpvf.p0_),
    n0_(tdpvf.n0_),
    U0_(tdpvf.U0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void planeContactDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void planeContactDisplacementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);

//     const planeContactDisplacementFvPatchVectorField& dmptf =
//         refCast
//         <
//             const planeContactDisplacementFvPatchVectorField
//         >(ptf);
}


// Update the coefficients associated with the patch field
void planeContactDisplacementFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    word UName = this->dimensionedInternalField().name();

    const fvsPatchField<tensor>& gradU =
        patch().lookupPatchField<surfaceTensorField, tensor>
        (
            "grad" + UName + "f"
        );

    const pointVectorField& pointU =
        mesh.lookupObject<pointVectorField>("point" + UName);

    const labelList& patchMeshPoints =
        mesh.boundaryMesh()[patch().index()].meshPoints();

    vectorField curLocalPoints =
        mesh.boundaryMesh()[patch().index()].localPoints();

    forAll(curLocalPoints, pointI)
    {
        curLocalPoints[pointI] += pointU[patchMeshPoints[pointI]];
    }

    PrimitivePatch<face, List, pointField> curPatch
    (
        mesh.boundaryMesh()[patch().index()].localFaces(),
        curLocalPoints
    );

    const vectorField& curCf = curPatch.faceCentres();


    // Set refValue and valueFraction

    const vectorField& Cf = patch().Cf();

    vector p = p0_ + U0_*n0_;

    vectorField n = patch().nf();

    vectorField delta(Cf.size(), vector::zero);
    symmTensorField valueFraction(Cf.size(), symmTensor::zero);

    forAll(curCf, faceI)
    {
        vector d = (p - Cf[faceI]);
        d = n0_*(n0_ & d);

        vector curD = (p - curCf[faceI]);
        curD = n0_*(n0_ & curD);

        if ((curD&n0_) >= -SMALL)
        {
            delta[faceI] = d;
//             valueFraction[faceI] = sqr(n[faceI]);
            valueFraction[faceI] = sqr(n0_);
        }
    }

    this->refValue() = delta;
    this->valueFraction() = valueFraction;

    // Set refGrad

    const fvPatchField<scalar>& mu =
        patch().lookupPatchField<volScalarField, scalar>("mu");

    const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>("lambda");

    Switch largeStrain(false);
//     (
//         patch().boundaryMesh().mesh().solutionDict().subDict("stressedFoam")
//        .lookup("largeStrain")
//     );

    vectorField traction(patch().size(), vector::zero);

    symmTensorField epsilon = symm(gradU);

    if (largeStrain)
    {
        epsilon += 0.5*symm(gradU & gradU.T());
    }

    symmTensorField sigma =
        2*mu*epsilon + I*(lambda*tr(epsilon));

    if (largeStrain)
    {
        tensorField F = I + gradU.T();

        tensorField invF = inv(F);

        scalarField SoS0 = mag(n & invF)*det(F);

        vectorField nCurrent = (n&invF);
        nCurrent /= mag(nCurrent);

        sigma = (1.0/det(F)) * symm(F & sigma & F.T());

        traction = nCurrent*(nCurrent&(nCurrent&sigma));

        traction = (invF & traction)*SoS0;
    }
    else
    {
//         traction = n*(n&(n&sigma));
        traction = -n0_*(n0_&(n0_&sigma));
    }

    forAll(traction, faceI)
    {
        if ( mag(this->valueFraction()[faceI]) < SMALL )
        {
            traction = vector::zero;
        }
    }

    refGrad() =
        traction
      - (n & (mu*gradU.T() - (mu + lambda)*gradU))
      - n*lambda*tr(gradU);

    if(largeStrain)
    {
        refGrad() -=
            (n & (mu*(gradU&gradU.T())))
          + 0.5*n*lambda*tr(gradU&gradU.T());
    }

    refGrad() /= (2.0*mu + lambda);














//     const fvPatchField<tensor>& gradDU =
//         patch().lookupPatchField<volTensorField, tensor>("grad(DU)");

//     tensorField DF = gradDU.T();

//     tensorField F = I + DF;

//     tensorField invF = inv(F);

//     vectorField n = patch().nf();

//     scalarField SoS0 = mag(n & invF)*det(F);

//     const fvPatchField<symmTensor>& sigma =
//         patch().lookupPatchField<volSymmTensorField, symmTensor>("sigma");

//     vectorField traction(Cf.size(), vector::zero);

//     if (max(mag(this->valueFraction())) > SMALL)
//     {
//         vectorField nCurrent = (n&invF);
//         nCurrent /= mag(nCurrent);

//         symmTensorField Depsilon = symm(gradDU)
//           + 0.5*symm(gradDU & gradDU.T());

//         symmTensorField DSigma =
//             2*mu*Depsilon + I*(lambda*tr(Depsilon));

//         symmTensorField newSigma = sigma + DSigma;

//         newSigma  = (1.0/det(F)) * symm(F & newSigma & F.T());

//         vectorField t = nCurrent*(nCurrent&(nCurrent&newSigma));

//         forAll(traction, faceI)
//         {
//             if (mag(this->valueFraction()[faceI]) > SMALL)
//             {
//                 traction[faceI] = t[faceI];
//             }
//         }
//     }

//     vectorField DTraction = (invF&traction)*SoS0 - (n & sigma);

//     vectorField newGradient =
//         DTraction
//       - (n & (mu*gradDU.T() - (mu + lambda)*gradDU))
//       - n*lambda*tr(gradDU);

//     newGradient -=
//         (n & (mu*(gradDU & gradDU.T())))
//       + 0.5*n*lambda*tr(gradDU & gradDU.T());

//     newGradient /= (2.0*mu + lambda);

//     refGrad() = newGradient;

    directionMixedFvPatchVectorField::updateCoeffs();
}


void planeContactDisplacementFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<vector> normalValue = transform(valueFraction(), refValue());

    word UName = this->dimensionedInternalField().name();

    //- non-ortho corrected gradValue
    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + UName + ")"
        );

    vectorField n = patch().nf();
    vectorField delta = patch().delta();
    vectorField k = delta - n*(n&delta);
    Field<vector> gradValue =
        this->patchInternalField()
      + (k&gradU.patchInternalField())
      + refGrad()/this->patch().deltaCoeffs();

    Field<vector> transformGradValue =
        transform(I - valueFraction(), gradValue);

    Field<vector>::operator=(normalValue + transformGradValue);

    fvPatchField<vector>::evaluate();
}

Foam::tmp<Foam::Field<vector> >
planeContactDisplacementFvPatchVectorField::snGrad() const
{
    Field<vector> normalValue = transform(valueFraction(), refValue());

    word UName = this->dimensionedInternalField().name();

    //- non-ortho corrected gradValue
    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + UName + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    //- correction vector
    vectorField k = delta - n*(n&delta);

    Field<vector> gradValue =
        this->patchInternalField()
      + (k&gradU.patchInternalField())
      + refGrad()/this->patch().deltaCoeffs();

    Field<vector> transformGradValue =
        transform(I - valueFraction(), gradValue);

    Field<vector> patchValue = normalValue + transformGradValue;

    // return
    // (
    //     *this
    //   - (patchInternalField() + (k&gradU.patchInternalField()))
    //   )*this->patch().deltaCoeffs();

    return
    (
        patchValue
      - (patchInternalField() + (k&gradU.patchInternalField()))
    )*this->patch().deltaCoeffs();
}


// Write
void planeContactDisplacementFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
//     traction_.writeEntry("traction", os);
//     pressure_.writeEntry("pressure", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, planeContactDisplacementFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
