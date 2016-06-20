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

#include "fixedHydPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "constitutiveModel.H"
#include "solidSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedHydPressureFvPatchScalarField::
fixedHydPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{
    fvPatchScalarField::operator==(patchInternalField());
}


fixedHydPressureFvPatchScalarField::
fixedHydPressureFvPatchScalarField
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


fixedHydPressureFvPatchScalarField::
fixedHydPressureFvPatchScalarField
(
    const fixedHydPressureFvPatchScalarField& tdpvf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(tdpvf, p, iF, mapper)
{}


fixedHydPressureFvPatchScalarField::
fixedHydPressureFvPatchScalarField
(
    const fixedHydPressureFvPatchScalarField& tdpvf
)
:
    fixedValueFvPatchScalarField(tdpvf)
{}


fixedHydPressureFvPatchScalarField::
fixedHydPressureFvPatchScalarField
(
    const fixedHydPressureFvPatchScalarField& tdpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tdpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedHydPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedHydPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void fixedHydPressureFvPatchScalarField::updateCoeffs()
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

    word DName = this->dimensionedInternalField().name();

    const fvsPatchField<tensor>& gradDf =
        patch().lookupPatchField<surfaceTensorField, tensor>
        (
            "gradDf"
        );

    const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>
        (
            "lambda"
        );

    // const fvPatchField<scalar>& mu =
    //     patch().lookupPatchField<volScalarField, scalar>
    //     (
    //         "mu"
    //     );

//     symmTensorField E = symm(gradDf);
//     if(nonLinear && !enforceLinear)
//     {
//         E += 0.5*symm(gradDf & gradDf.T());
//     }

//     fvPatchField<scalar>::operator==(-(lambda+(2.0/3.0)*mu)*tr(E));

    scalarField Jf = det(I + gradDf.T());
    scalarField e0f = (Jf - 1.0)/3.0;
    fvPatchField<scalar>::operator==(-3*lambda*(e0f + 0.5*sqr(e0f)));

//     Info << "Update hyd pressure at " << patch().name() << endl;

    fixedValueFvPatchScalarField::updateCoeffs();
}


Foam::tmp<Foam::Field<scalar> > fixedHydPressureFvPatchScalarField::
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

tmp<Field<scalar> > fixedHydPressureFvPatchScalarField::
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
void fixedHydPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, fixedHydPressureFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
