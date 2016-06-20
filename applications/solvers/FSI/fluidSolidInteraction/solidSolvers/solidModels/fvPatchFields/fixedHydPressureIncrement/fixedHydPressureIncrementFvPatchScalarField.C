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

#include "fixedHydPressureIncrementFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "constitutiveModel.H"
#include "solidSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedHydPressureIncrementFvPatchScalarField::
fixedHydPressureIncrementFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{
    fvPatchScalarField::operator==(patchInternalField());
}


fixedHydPressureIncrementFvPatchScalarField::
fixedHydPressureIncrementFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF)
{
    fvPatchScalarField::operator==(patchInternalField());

    Info << "Creating fixed hyd pressure boundary conditions" << endl;
}


fixedHydPressureIncrementFvPatchScalarField::
fixedHydPressureIncrementFvPatchScalarField
(
    const fixedHydPressureIncrementFvPatchScalarField& tdpvf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(tdpvf, p, iF, mapper)
{}


fixedHydPressureIncrementFvPatchScalarField::
fixedHydPressureIncrementFvPatchScalarField
(
    const fixedHydPressureIncrementFvPatchScalarField& tdpvf
)
:
    fixedValueFvPatchScalarField(tdpvf)
{}


fixedHydPressureIncrementFvPatchScalarField::
fixedHydPressureIncrementFvPatchScalarField
(
    const fixedHydPressureIncrementFvPatchScalarField& tdpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tdpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedHydPressureIncrementFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedHydPressureIncrementFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void fixedHydPressureIncrementFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

//     // Looking up solid solver
//     const solidSolver& stress =
//         this->db().objectRegistry::lookupObject<solidSolver>
//         (
//             "solidProperties"
//         );

//     Switch nonLinear
//     (
//         stress.solidProperties().lookup("nonLinear")
//     );

//     Switch enforceLinear
//     (
//         stress.solidProperties().lookup("enforceLinear")
//     );

//     word DName = this->dimensionedInternalField().name();

//     const fvsPatchField<tensor>& gradDf =
//         patch().lookupPatchField<surfaceTensorField, tensor>
//         (
//             "gradDf"
//         );

    const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>
        (
            "lambda"
        );

    const fvPatchField<scalar>& mu =
        patch().lookupPatchField<volScalarField, scalar>
        (
            "mu"
        );

    scalarField K = (lambda + (2.0/3.0)*mu);

    const fvsPatchField<scalar>& J =
        patch().lookupPatchField<surfaceScalarField, tensor>
        (
            "Jf"
        );


//     symmTensorField E = symm(gradDf);
//     if(nonLinear && !enforceLinear)
//     {
//         E += 0.5*symm(gradDf & gradDf.T());
//     }

//     fvPatchField<scalar>::operator==(-(lambda+(2.0/3.0)*mu)*tr(E));

//     scalarField Jf = det(I + gradDf.T());
//     scalarField e0f = (Jf - 1.0)/3.0;

    fvPatchField<scalar>::operator==(0.5*K*(sqr(J)-1));

//     Info << "Update hyd pressure at " << patch().name() << endl;

    fixedValueFvPatchScalarField::updateCoeffs();
}


Foam::tmp<Foam::Field<scalar> > fixedHydPressureIncrementFvPatchScalarField::
snGrad() const
{
    word pName = this->dimensionedInternalField().name();

    const fvPatchField<vector>& gradp =
        patch().lookupPatchField<volVectorField, vector>
        (
            "grad(" + pName + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

//     return
//     (
//         *this
//       - (this->patchInternalField() + (k&gradp.patchInternalField()))
//     )*this->patch().deltaCoeffs();

    scalarField dpP = (k&gradp.patchInternalField());
    scalarField nGradpP = (n&gradp.patchInternalField());

    return
        2
       *(
            *this
          - (this->patchInternalField() + dpP)
        )*this->patch().deltaCoeffs()
      - nGradpP;
}

tmp<Field<scalar> > fixedHydPressureIncrementFvPatchScalarField::
gradientBoundaryCoeffs() const
{
    word pName = this->dimensionedInternalField().name();

    const fvPatchField<vector>& gradp =
        patch().lookupPatchField<volVectorField, vector>
        (
            "grad(" + pName + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    vectorField k = delta - n*(n&delta);

//     return this->patch().deltaCoeffs()
//        *(*this - (k&gradp.patchInternalField()));

    scalarField dpP = (k&gradp.patchInternalField());
    scalarField nGradpP = (n&gradp.patchInternalField());

    return
        this->patch().deltaCoeffs()
       *(
           2*(*this - dpP)
         - this->patchInternalField()
        )
      - nGradpP;
}


// Write
void fixedHydPressureIncrementFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    fixedHydPressureIncrementFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
