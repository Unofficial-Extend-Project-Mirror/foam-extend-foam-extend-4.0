/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "solidTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "tractionBoundaryGradient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    fieldName_("undefined"),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    nonLinear_(nonLinearGeometry::OFF),
    orthotropic_(false)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    fieldName_(dimensionedInternalField().name()),
    traction_("traction", dict, p.size()),
    pressure_("pressure", dict, p.size()),
    nonLinear_(nonLinearGeometry::OFF),
    orthotropic_(false)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;

    Info << "Patch " << patch().name()
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
            Info << "\tnonLinear set to updated Lagrangian"
                << endl;
        }
        else if (nonLinear_ == nonLinearGeometry::TOTAL_LAGRANGIAN)
        {
            Info << "\tnonLinear set to total Lagrangian"
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


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(stpvf, p, iF, mapper),
    fieldName_(stpvf.fieldName_),
    traction_(stpvf.traction_, mapper),
    pressure_(stpvf.pressure_, mapper),
    nonLinear_(stpvf.nonLinear_),
    orthotropic_(stpvf.orthotropic_)
{}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& stpvf
)
:
    fixedGradientFvPatchVectorField(stpvf),
    fieldName_(stpvf.fieldName_),
    traction_(stpvf.traction_),
    pressure_(stpvf.pressure_),
    nonLinear_(stpvf.nonLinear_),
    orthotropic_(stpvf.orthotropic_)
{}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(stpvf, iF),
    fieldName_(stpvf.fieldName_),
    traction_(stpvf.traction_),
    pressure_(stpvf.pressure_),
    nonLinear_(stpvf.nonLinear_),
    orthotropic_(stpvf.orthotropic_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void solidTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const solidTractionFvPatchVectorField& dmptf =
        refCast<const solidTractionFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void solidTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    gradient() = tractionBoundaryGradient()
      (
       traction_,
       pressure_,
       word(fieldName_),
       patch(),
       orthotropic_,
       nonLinearGeometry::nonLinearNames_[nonLinear_]
       )();

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void solidTractionFvPatchVectorField::evaluate(const Pstream::commsTypes)
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

    //- non-orthogonal correction vectors
    vectorField k = delta - n*(n&delta);

    Field<vector>::operator=
    (
        this->patchInternalField()
      + (k&gradField.patchInternalField())
      + gradient()/this->patch().deltaCoeffs()
    );

    fvPatchField<vector>::evaluate();
}

// Write
void solidTractionFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("nonLinear")
        << nonLinearGeometry::nonLinearNames_[nonLinear_]
        << token::END_STATEMENT << nl;

    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidTractionFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
