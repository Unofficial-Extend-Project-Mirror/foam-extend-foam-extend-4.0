/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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

Class
    timeVaryingFixedDisplacementZeroShearFvPatchVectorField

Description

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

timeVaryingFixedDisplacementZeroShearFvPatchVectorField::timeVaryingFixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    fieldName_("undefined"),
    nonLinear_(OFF),
    orthotropic_(false),
    timeSeries_()
{}


timeVaryingFixedDisplacementZeroShearFvPatchVectorField::timeVaryingFixedDisplacementZeroShearFvPatchVectorField
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


timeVaryingFixedDisplacementZeroShearFvPatchVectorField::timeVaryingFixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    fieldName_(dimensionedInternalField().name()),
    nonLinear_(OFF),
    orthotropic_(false),
    timeSeries_(dict)
{
    //- check if traction boundary is for non linear solver
    if(dict.found("nonLinear"))
      {
	nonLinear_ = nonLinearNames_.read(dict.lookup("nonLinear"));;

	if(nonLinear_ == UPDATED_LAGRANGIAN)
	  {
	    Info << "\tnonLinear set to updated Lagrangian"
		 << endl;
	  }
	else if(nonLinear_ == TOTAL_LAGRANGIAN)
	  {
	    Info << "\tnonLinear set to total Lagrangian"
		 << endl;
	  }
      }

    if(dict.found("orthotropic"))
      {
	orthotropic_ = Switch(dict.lookup("orthotropic"));
	Info << "\t\torthotropic set to " << orthotropic_ << endl;
      }

    //- the leastSquares has zero non-orthogonal correction
    //- on the boundary
    //- so the gradient scheme should be extendedLeastSquares
//     if(Foam::word(dimensionedInternalField().mesh().gradScheme("grad(" + fieldName_ + ")")) != "extendedLeastSquares")
//       {
// 	Warning << "The gradScheme for " << fieldName_
// 		<< " should be \"extendedLeastSquares 0\" for the boundary "
// 		<< "non-orthogonal correction to be right" << endl;
//       }

  this->refGrad() = vector::zero;
  
  vectorField n = patch().nf();      
  this->valueFraction() = sqr(n);

  if (dict.found("value"))
    {
      Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
  else
    {
      FatalError << "value entry not found for patch " << patch().name() << endl;
    }
  //this->refValue() = *this;
  this->refValue() = timeSeries_(this->db().time().timeOutputValue());

  Field<vector> normalValue = transform(valueFraction(), refValue());
  
  Field<vector> gradValue =
    this->patchInternalField() + refGrad()/this->patch().deltaCoeffs();
  
  Field<vector> transformGradValue =
    transform(I - valueFraction(), gradValue);
  
  Field<vector>::operator=(normalValue + transformGradValue);
}


timeVaryingFixedDisplacementZeroShearFvPatchVectorField::timeVaryingFixedDisplacementZeroShearFvPatchVectorField
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
    vectorField disp(patch().size(), timeSeries_(this->db().time().timeOutputValue()));
    if(fieldName_ == "DU")
      {
        const fvPatchField<vector>& U =
          patch().lookupPatchField<volVectorField, vector>("U");
	disp -= U;
      }
    else if(fieldName_ != "U")
      {
        FatalError << "The displacement field should be U or DU"
                   << exit(FatalError);
      }
    this->refValue() = disp;


    // set value fraction to fix reference patch normal
    // only done at initialisation above
    //vectorField n = patch().nf();
    //this->valueFraction() = sqr(n);

    refGrad() = tractionBoundaryGradient()
      (
       vectorField(patch().size(), vector::zero),
       scalarField(patch().size(), 0.0),
       word(fieldName_),
       patch(),
       orthotropic_,
       NamedEnum<Foam::timeVaryingFixedDisplacementZeroShearFvPatchVectorField::nonLinearType, 3>::names[nonLinear_]
       )();

    directionMixedFvPatchVectorField::updateCoeffs();
}


// Write
void timeVaryingFixedDisplacementZeroShearFvPatchVectorField::write(Ostream& os) const
{
    directionMixedFvPatchVectorField::write(os);
    os.writeKeyword("nonLinear") << nonLinearNames_[nonLinear_] << token::END_STATEMENT << nl;
    timeSeries_.write(os);
}


template<>
const char* Foam::NamedEnum<Foam::timeVaryingFixedDisplacementZeroShearFvPatchVectorField::nonLinearType, 3>::names[] =
  {
    "off",
    "updatedLagrangian",
    "totalLagrangian"
  };

const Foam::NamedEnum<Foam::timeVaryingFixedDisplacementZeroShearFvPatchVectorField::nonLinearType, 3>
Foam::timeVaryingFixedDisplacementZeroShearFvPatchVectorField::nonLinearNames_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, timeVaryingFixedDisplacementZeroShearFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
