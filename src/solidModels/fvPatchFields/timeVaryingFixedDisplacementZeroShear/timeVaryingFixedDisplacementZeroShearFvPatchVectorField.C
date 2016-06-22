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

Class
    timeVaryingFixedDisplacementZeroShearFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "timeVaryingFixedDisplacementZeroShearFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "tractionBoundaryGradient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

timeVaryingFixedDisplacementZeroShearFvPatchVectorField::
timeVaryingFixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    fieldName_("undefined"),
    nonLinear_(nonLinearGeometry::OFF),
    orthotropic_(false),
    timeSeries_()
{}


timeVaryingFixedDisplacementZeroShearFvPatchVectorField::
timeVaryingFixedDisplacementZeroShearFvPatchVectorField
(
    const timeVaryingFixedDisplacementZeroShearFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_),
    nonLinear_(ptf.nonLinear_),
    orthotropic_(ptf.orthotropic_),
    timeSeries_(ptf.timeSeries_)
{}


timeVaryingFixedDisplacementZeroShearFvPatchVectorField::
timeVaryingFixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    fieldName_(dimensionedInternalField().name()),
    nonLinear_(nonLinearGeometry::OFF),
    orthotropic_(false),
    timeSeries_(dict)
{
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
        Info << "\t\torthotropic set to " << orthotropic_ << endl;
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

    this->refGrad() = vector::zero;

    vectorField n = patch().nf();
    this->valueFraction() = sqr(n);

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        FatalError
            << "value entry not found for patch " << patch().name() << endl;
    }

    this->refValue() = timeSeries_(this->db().time().timeOutputValue());

    Field<vector> normalValue = transform(valueFraction(), refValue());

    Field<vector> gradValue =
        this->patchInternalField() + refGrad()/this->patch().deltaCoeffs();

    Field<vector> transformGradValue =
        transform(I - valueFraction(), gradValue);

    Field<vector>::operator=(normalValue + transformGradValue);
}


timeVaryingFixedDisplacementZeroShearFvPatchVectorField::
timeVaryingFixedDisplacementZeroShearFvPatchVectorField
(
    const timeVaryingFixedDisplacementZeroShearFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF),
    fieldName_(ptf.fieldName_),
    nonLinear_(ptf.nonLinear_),
    orthotropic_(ptf.orthotropic_),
    timeSeries_(ptf.timeSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void timeVaryingFixedDisplacementZeroShearFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void timeVaryingFixedDisplacementZeroShearFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);
}


void timeVaryingFixedDisplacementZeroShearFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // set refValue
    vectorField disp
    (
        patch().size(),
        timeSeries_(this->db().time().timeOutputValue())
    );

    if (fieldName_ == "DU")
    {
          const fvPatchField<vector>& U =
              patch().lookupPatchField<volVectorField, vector>("U");
          disp -= U;
    }
    else if (fieldName_ != "U")
    {
        FatalError
            << "The displacement field should be U or DU"
            << exit(FatalError);
    }
    this->refValue() = disp;

    // set value fraction to fix reference patch normal
    // only done at initialisation above
    //vectorField n = patch().nf();
    //this->valueFraction() = sqr(n);

    bool incremental(fieldName_ == "DU");

    refGrad() = tractionBoundaryGradient::snGrad
    (
        vectorField(patch().size(), vector::zero),
        scalarField(patch().size(), 0),
        fieldName_,
        "U",
        patch(),
        orthotropic_,
        nonLinear_,
        incremental
    );

    directionMixedFvPatchVectorField::updateCoeffs();
}


// Write
void
timeVaryingFixedDisplacementZeroShearFvPatchVectorField::write
(
    Ostream& os
) const
{
    directionMixedFvPatchVectorField::write(os);
    os.writeKeyword("nonLinear")
        << nonLinearGeometry::nonLinearNames_[nonLinear_]
        << token::END_STATEMENT << nl;
    timeSeries_.write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    timeVaryingFixedDisplacementZeroShearFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
