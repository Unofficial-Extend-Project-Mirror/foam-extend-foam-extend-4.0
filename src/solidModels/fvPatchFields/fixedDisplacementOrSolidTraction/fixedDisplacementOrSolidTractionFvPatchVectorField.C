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
    fixedDisplacementOrSolidTractionFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "fixedDisplacementOrSolidTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "tractionBoundaryGradient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedDisplacementOrSolidTractionFvPatchVectorField::
fixedDisplacementOrSolidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    fieldName_("undefined"),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    displacement_(p.size(), vector::zero),
    nonLinear_(OFF),
    orthotropic_(false),
    timeSeries_()
{}


fixedDisplacementOrSolidTractionFvPatchVectorField::
fixedDisplacementOrSolidTractionFvPatchVectorField
(
    const fixedDisplacementOrSolidTractionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_),
    traction_(ptf.traction_, mapper),
    pressure_(ptf.pressure_, mapper),
    displacement_(ptf.displacement_, mapper),
    nonLinear_(ptf.nonLinear_),
    orthotropic_(ptf.orthotropic_),
    timeSeries_(ptf.timeSeries_)
{}


fixedDisplacementOrSolidTractionFvPatchVectorField::
fixedDisplacementOrSolidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    fieldName_(dimensionedInternalField().name()),
    traction_("traction", dict, p.size()),
    pressure_("pressure", dict, p.size()),
    displacement_("displacement", dict, p.size()),
    nonLinear_(OFF),
    orthotropic_(false),
    timeSeries_(dict)
{
  Info<< "fixedDisplacementOrSolidTraction boundary condition"
       << endl;

    //- check if traction boundary is for non linear solver
    if (dict.found("nonLinear"))
      {
          nonLinear_ = nonLinearNames_.read(dict.lookup("nonLinear"));

          if (nonLinear_ == UPDATED_LAGRANGIAN)
          {
              Info<< "\tnonLinear set to updated Lagrangian"
                   << endl;
          }
          else if (nonLinear_ == TOTAL_LAGRANGIAN)
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

  this->refValue() = displacement_;
  this->refGrad() = vector::zero;
  this->valueFraction() = symmTensor(1,0,0,1,0,1);

  Field<vector>::operator=(displacement_);
}


fixedDisplacementOrSolidTractionFvPatchVectorField::
fixedDisplacementOrSolidTractionFvPatchVectorField
(
    const fixedDisplacementOrSolidTractionFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF),
    fieldName_(ptf.fieldName_),
    traction_(ptf.traction_),
    pressure_(ptf.pressure_),
    displacement_(ptf.displacement_),
    nonLinear_(ptf.nonLinear_),
    orthotropic_(ptf.orthotropic_),
    timeSeries_(ptf.timeSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void fixedDisplacementOrSolidTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedDisplacementOrSolidTractionFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);
}


void fixedDisplacementOrSolidTractionFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if ( mag(timeSeries_(this->db().time().timeOutputValue())) < SMALL)
    {
        // traction boundary

        // set valueFraction to zero
        this->valueFraction() = symmTensor::zero;

        // set gradient to enfore specified traction
        refGrad() = tractionBoundaryGradient::snGrad
        (
            traction_,
            pressure_,
            fieldName_,
            "U",
            patch(),
            orthotropic_,
            nonLinearGeometry::nonLinearNames_[nonLinear_]
        );
    }
    else
    {
        // fixed displacement

        // set valueFraction to one
        this->valueFraction() = symmTensor(1,0,0,1,0,1);

        // set displacement
        refValue() = displacement_;
      }

    directionMixedFvPatchVectorField::updateCoeffs();
}


// Write
void fixedDisplacementOrSolidTractionFvPatchVectorField::write
(Ostream& os) const
{
    directionMixedFvPatchVectorField::write(os);
    os.writeKeyword("nonLinear")
        << nonLinearGeometry::nonLinearNames_[nonLinear_]
        << token::END_STATEMENT << nl;
    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);
    displacement_.writeEntry("displacement", os);
    timeSeries_.write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedDisplacementOrSolidTractionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
