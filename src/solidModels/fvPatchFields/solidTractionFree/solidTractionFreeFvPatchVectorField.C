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

#include "solidTractionFreeFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "tractionBoundaryGradient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidTractionFreeFvPatchVectorField::
solidTractionFreeFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    fieldName_("undefined"),
    nonLinear_(nonLinearGeometry::OFF),
    orthotropic_(false)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


solidTractionFreeFvPatchVectorField::
solidTractionFreeFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    fieldName_(dimensionedInternalField().name()),
    nonLinear_(nonLinearGeometry::OFF),
    orthotropic_(false)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;

    Info<< "Patch " << patch().name()
        << "\tTraction boundary field: " << fieldName_ << endl;

    //- check if traction boundary is for non linear solver
    if (dict.found("nonLinear"))
    {
        nonLinear_ = nonLinearGeometry::nonLinearNames_.read
        (
            dict.lookup("nonLinear")
        );

        if (nonLinear_ == nonLinearGeometry::UPDATED_LAGRANGIAN)
        {
            Info<< "\tnonLinear set to updated Lagrangian"
                << endl;
        }
        else if (nonLinear_ == nonLinearGeometry::TOTAL_LAGRANGIAN)
        {
            Info<< "\tnonLinear set to total Lagrangian"
                << endl;
        }
    }

    if (dict.found("orthotropic"))
    {
        orthotropic_ = Switch(dict.lookup("orthotropic"));
        Info<< "\t\torthotropic set to " << orthotropic_ << endl;
    }

    //- the leastSquares has zero non-orthogonal correction
    //- on the boundary
    //- so the gradient scheme should be extendedLeastSquares
    if
    (
        Foam::word
        (
            dimensionedInternalField().mesh().schemesDict().gradScheme
            (
                "grad(" + fieldName_ + ")"
            )
        ) != "extendedLeastSquares"
    )
    {
        Warning << "The gradScheme for " << fieldName_
            << " should be \"extendedLeastSquares 0\" for the boundary "
            << "non-orthogonal correction to be right" << endl;
    }
}


solidTractionFreeFvPatchVectorField::
solidTractionFreeFvPatchVectorField
(
    const solidTractionFreeFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(stpvf, p, iF, mapper),
    fieldName_(stpvf.fieldName_),
    nonLinear_(stpvf.nonLinear_),
    orthotropic_(stpvf.orthotropic_)
{}


solidTractionFreeFvPatchVectorField::
solidTractionFreeFvPatchVectorField
(
    const solidTractionFreeFvPatchVectorField& stpvf
)
:
    fixedGradientFvPatchVectorField(stpvf),
    fieldName_(stpvf.fieldName_),
    nonLinear_(stpvf.nonLinear_),
    orthotropic_(stpvf.orthotropic_)
{}


solidTractionFreeFvPatchVectorField::
solidTractionFreeFvPatchVectorField
(
    const solidTractionFreeFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(stpvf, iF),
    fieldName_(stpvf.fieldName_),
    nonLinear_(stpvf.nonLinear_),
    orthotropic_(stpvf.orthotropic_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidTractionFreeFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void solidTractionFreeFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
void solidTractionFreeFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    bool incremental(fieldName_ == "DU");

    gradient() = tractionBoundaryGradient::snGrad
    (
        vectorField(patch().size(), vector::zero),
        scalarField(patch().size(), 0.0),
        fieldName_,
        "U",
        patch(),
        orthotropic_,
        nonLinear_,
        incremental
    );

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void solidTractionFreeFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + fieldName_ + ")"
        );

    vectorField n = patch().nf();
    vectorField delta = patch().delta();

    vectorField k = delta - n*(n&delta);

    Field<vector>::operator=
    (
        this->patchInternalField()
      + (k & gradField.patchInternalField())
      + gradient()/this->patch().deltaCoeffs()
    );

    fvPatchField<vector>::evaluate();
}

// Write
void solidTractionFreeFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("nonLinear")
        << nonLinearGeometry::nonLinearNames_[nonLinear_]
    << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidTractionFreeFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
