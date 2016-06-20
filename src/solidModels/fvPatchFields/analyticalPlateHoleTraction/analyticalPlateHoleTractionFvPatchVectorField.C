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

#include "analyticalPlateHoleTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "constitutiveModel.H"
#include "volFields.H"
#include "fvc.H"
#include "fixedValueFvPatchFields.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

analyticalPlateHoleTractionFvPatchVectorField::
analyticalPlateHoleTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    UName_("undefined")
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


analyticalPlateHoleTractionFvPatchVectorField::
analyticalPlateHoleTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    UName_("U")
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


analyticalPlateHoleTractionFvPatchVectorField::
analyticalPlateHoleTractionFvPatchVectorField
(
    const analyticalPlateHoleTractionFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(stpvf, p, iF, mapper),
    UName_(stpvf.UName_)
{}


analyticalPlateHoleTractionFvPatchVectorField::
analyticalPlateHoleTractionFvPatchVectorField
(
    const analyticalPlateHoleTractionFvPatchVectorField& stpvf
)
:
    fixedGradientFvPatchVectorField(stpvf),
    UName_(stpvf.UName_)
{}


analyticalPlateHoleTractionFvPatchVectorField::
analyticalPlateHoleTractionFvPatchVectorField
(
    const analyticalPlateHoleTractionFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(stpvf, iF),
    UName_(stpvf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void analyticalPlateHoleTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void analyticalPlateHoleTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void analyticalPlateHoleTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    vectorField n = patch().nf();

     const constitutiveModel& rheology =
         this->db().objectRegistry::lookupObject<constitutiveModel>
         ("rheologyProperties");
    scalarField mu =
        rheology.mu()().boundaryField()[patch().index()];
    scalarField lambda =
        rheology.lambda()().boundaryField()[patch().index()];

    const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>
        ("grad(" + UName_ + ")");

    vectorField Traction(n.size(),vector::zero);

    const vectorField& Cf = patch().Cf();

    forAll(Traction, faceI)
      {
    vector curC(Cf[faceI].x(), Cf[faceI].y(), 0);
    vector curN = n[faceI];

    if (patch().name() == "hole")
      {
        curC /= mag(curC);
        curC *= 0.5;

        curN = -curC/mag(curC);
      }

    Traction[faceI] =
      (n[faceI] & plateHoleSolution(curC));
      }

    //- set patch gradient
    vectorField newGradient =
      Traction
      - (n & (mu*gradU.T() - (mu + lambda)*gradU))
      - n*lambda*tr(gradU);

    newGradient /= (2.0*mu + lambda);

    gradient() = newGradient;

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void analyticalPlateHoleTractionFvPatchVectorField::evaluate
(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const fvPatchField<tensor>& gradDU =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + UName_ + ")"
        );

    vectorField n = patch().nf();
    vectorField delta = patch().delta();

    vectorField k = delta - n*(n&delta);

    Field<vector>::operator=
    (
        this->patchInternalField()
      + (k&gradDU.patchInternalField())
      + gradient()/this->patch().deltaCoeffs()
    );

    fvPatchField<vector>::evaluate();
}

// Write
void analyticalPlateHoleTractionFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}

symmTensor analyticalPlateHoleTractionFvPatchVectorField::plateHoleSolution
(const vector& C)
{
    tensor sigma = tensor::zero;

    scalar T = 10000;
    scalar a = 0.5;

    scalar r = ::sqrt(sqr(C.x()) + sqr(C.y()));
    scalar theta = Foam::atan2(C.y(), C.x());

    coordinateSystem cs("polarCS", C, vector(0, 0, 1), C/mag(C));

    sigma.xx() =
        T*(1 - sqr(a)/sqr(r))/2
      + T*(1 + 3*pow(a,4)/pow(r,4) - 4*sqr(a)/sqr(r))*::cos(2*theta)/2;

    sigma.xy() =
      - T*(1 - 3*pow(a,4)/pow(r,4) + 2*sqr(a)/sqr(r))*::sin(2*theta)/2;

    sigma.yx() = sigma.xy();

    sigma.yy() =
        T*(1 + sqr(a)/sqr(r))/2
      - T*(1 + 3*pow(a,4)/pow(r,4))*::cos(2*theta)/2;


    // Transformation to global coordinate system
    sigma = ((cs.R()&sigma)&cs.R().T());

    symmTensor S = symmTensor::zero;

    S.xx() = sigma.xx();
    S.xy() = sigma.xy();
    S.yy() = sigma.yy();

    return S;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    analyticalPlateHoleTractionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
