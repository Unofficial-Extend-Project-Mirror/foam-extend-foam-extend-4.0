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

#include "directionMixedDisplacementFvPatchVectorField.H"
#include "symmTransformField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "solidSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

directionMixedDisplacementFvPatchVectorField
::directionMixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    limitCoeff_(1.0)
{
    // Looking up solid solver
    if
    (
        this->db().objectRegistry::foundObject<solidSolver>
        (
            "solidProperties"
        )
    )
    {
        const solidSolver& stress =
            this->db().objectRegistry::lookupObject<solidSolver>
            (
                "solidProperties"
            );

        if (stress.solidProperties().found("snGradLimitCoeff"))
        {
            limitCoeff_ =
                readScalar
                (
                    stress.solidProperties().lookup("snGradLimitCoeff")
                );

            Info << "snGradLimitCoeff: " << limitCoeff_ << endl;
        }
    }
}


directionMixedDisplacementFvPatchVectorField
::directionMixedDisplacementFvPatchVectorField
(
    const directionMixedDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    limitCoeff_(ptf.limitCoeff_)
{}


directionMixedDisplacementFvPatchVectorField
::directionMixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF, dict),
    limitCoeff_(1.0)
{
    Info << "Direction mixed boundary condition with non-orthogonal correction"
        << endl;
    directionMixedFvPatchVectorField::evaluate();

    if (dict.found("limitCoeff"))
    {
        limitCoeff_ =
            scalar(readScalar(dict.lookup("limitCoeff")));
        Info << "Limiter coefficient: " << limitCoeff_ << endl;
    }

    // Looking up solid solver
    if
    (
        this->db().objectRegistry::foundObject<solidSolver>
        (
            "solidProperties"
        )
    )
    {
        const solidSolver& stress =
            this->db().objectRegistry::lookupObject<solidSolver>
            (
                "solidProperties"
            );

        if (stress.solidProperties().found("snGradLimitCoeff"))
        {
            limitCoeff_ =
                readScalar
                (
                    stress.solidProperties().lookup("snGradLimitCoeff")
                );

            Info << "snGradLimitCoeff: " << limitCoeff_ << endl;
        }
    }
}


directionMixedDisplacementFvPatchVectorField
::directionMixedDisplacementFvPatchVectorField
(
    const directionMixedDisplacementFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF),
    limitCoeff_(ptf.limitCoeff_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void directionMixedDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedFvPatchVectorField::autoMap(m);
}


void directionMixedDisplacementFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);

//     const directionMixedDisplacementFvPatchVectorField& dmptf =
//         refCast<const directionMixedDisplacementFvPatchVectorField >(ptf);
}


tmp<Field<vector> > directionMixedDisplacementFvPatchVectorField
::snGrad() const
{
    bool secondOrder_(false);

    Field<vector> pif = this->patchInternalField();

    Field<vector> normalValue =
        transform(this->valueFraction(), this->refValue());

    vectorField n = this->patch().nf();
    vectorField delta = patch().delta();
    vectorField k = delta - n*(n&delta);

    const fvPatchField<tensor>& gradD =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + this->dimensionedInternalField().name() + ")"
        );

    vectorField snGradCorrection =
      - (k&gradD.patchInternalField())
       *this->patch().deltaCoeffs();

    if (limitCoeff_ < 1-SMALL)
    {
        vectorField uncorrectedSnGrad =
            (
               *this
              - this->patchInternalField()
            )
           *this->patch().deltaCoeffs();

        scalarField limiter =
            (
                min
                (
                    limitCoeff_*mag(uncorrectedSnGrad + snGradCorrection)
                   /((1 - limitCoeff_)*mag(snGradCorrection) + SMALL),
                    1.0
                )
            );

        snGradCorrection *= limiter;
    }


    //
    Field<vector> gradValue =
        pif
      - snGradCorrection/this->patch().deltaCoeffs()
//       + (k&gradD.patchInternalField())
      + this->refGrad()/this->patch().deltaCoeffs();

    if (secondOrder_)
    {
        vectorField nGradDP = (n&gradD.patchInternalField());

        gradValue =
            this->patchInternalField() + (k&gradD.patchInternalField())
          + 0.5*(nGradDP + this->refGrad())/this->patch().deltaCoeffs();
    }

    Field<vector> transformGradValue =
        transform(I - this->valueFraction(), gradValue);

    if (secondOrder_)
    {
        vectorField nGradDP = (n&gradD.patchInternalField());

        return
            2
           *(
               normalValue + transformGradValue
             - (pif + (k&gradD.patchInternalField()))
           )*this->patch().deltaCoeffs()
         - nGradDP;
    }

    return
        (
            normalValue + transformGradValue
          - (
              pif
            - snGradCorrection/this->patch().deltaCoeffs()
//           + (k&gradD.patchInternalField())
            )
        )
       *this->patch().deltaCoeffs();
}


void directionMixedDisplacementFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    bool secondOrder_(false);

    vectorField n = this->patch().nf();
    vectorField delta = patch().delta();
    vectorField k = delta - n*(n&delta);

    const fvPatchField<tensor>& gradD =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + this->dimensionedInternalField().name() + ")"
        );


    // Calc limited snGrad correction

    vectorField snGradCorrection =
      - (k&gradD.patchInternalField())
       *this->patch().deltaCoeffs();

    if (limitCoeff_ < 1-SMALL)
    {
        vectorField uncorrectedSnGrad =
            (
               *this
              - this->patchInternalField()
            )
           *this->patch().deltaCoeffs();

        scalarField limiter =
            (
                min
                (
                    limitCoeff_*mag(uncorrectedSnGrad + snGradCorrection)
                   /((1 - limitCoeff_)*mag(snGradCorrection) + SMALL),
                    1.0
                )
            );

        snGradCorrection *= limiter;
    }

    //
    Field<vector> normalValue =
        transform(this->valueFraction(), this->refValue());

    Field<vector> gradValue =
        this->patchInternalField()
      - snGradCorrection/this->patch().deltaCoeffs()
//       + (k&gradD.patchInternalField())
      + this->refGrad()/this->patch().deltaCoeffs();

    if (secondOrder_)
    {
        vectorField nGradDP = (n&gradD.patchInternalField());

        gradValue =
            this->patchInternalField()
          + (k&gradD.patchInternalField())
          + 0.5*(nGradDP + this->refGrad())/this->patch().deltaCoeffs();
    }

    Field<vector> transformGradValue =
        transform(I - this->valueFraction(), gradValue);

    Field<vector>::operator=(normalValue + transformGradValue);

    fvPatchField<vector>::evaluate();
}


void directionMixedDisplacementFvPatchVectorField::write(Ostream& os) const
{
    directionMixedFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    directionMixedDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
