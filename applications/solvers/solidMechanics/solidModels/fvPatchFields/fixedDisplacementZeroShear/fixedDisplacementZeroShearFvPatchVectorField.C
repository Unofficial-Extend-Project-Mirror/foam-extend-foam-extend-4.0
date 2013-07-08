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
    fixedDisplacementZeroShearFvPatchVectorField

Description

\*---------------------------------------------------------------------------*/

#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "rheologyModel.H"
#include "plasticityModel.H"
#include "thermalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedDisplacementZeroShearFvPatchVectorField::fixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    fieldName_("undefined"),
    nonLinear_(OFF)
{}


fixedDisplacementZeroShearFvPatchVectorField::fixedDisplacementZeroShearFvPatchVectorField
(
    const fixedDisplacementZeroShearFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_),
    nonLinear_(ptf.nonLinear_)
{}


fixedDisplacementZeroShearFvPatchVectorField::fixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    fieldName_(dimensionedInternalField().name()),
    nonLinear_(OFF)
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

    //- the leastSquares has zero non-orthogonal correction
    //- on the boundary
    //- so the gradient scheme should be extendedLeastSquares
    if(Foam::word(dimensionedInternalField().mesh().gradScheme("grad(" + fieldName_ + ")")) != "extendedLeastSquares")
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
      FatalError << "value entry not found for patch " << patch().name() << endl;
    }
  this->refValue() = *this;

  Field<vector> normalValue = transform(valueFraction(), refValue());
  
  Field<vector> gradValue =
    this->patchInternalField() + refGrad()/this->patch().deltaCoeffs();
  
  Field<vector> transformGradValue =
    transform(I - valueFraction(), gradValue);
  
  Field<vector>::operator=(normalValue + transformGradValue);
}


fixedDisplacementZeroShearFvPatchVectorField::fixedDisplacementZeroShearFvPatchVectorField
(
    const fixedDisplacementZeroShearFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF),
    fieldName_(ptf.fieldName_),
    nonLinear_(ptf.nonLinear_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void fixedDisplacementZeroShearFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedDisplacementZeroShearFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);
}


void fixedDisplacementZeroShearFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //---------------------------//
    //- material properties
    //---------------------------//
    const rheologyModel& rheology =
        this->db().objectRegistry::lookupObject<rheologyModel>("rheologyProperties");
    scalarField mu = 
        rheology.mu()().boundaryField()[patch().index()];
    scalarField lambda =
        rheology.lambda()().boundaryField()[patch().index()];

    if(rheology.type() == plasticityModel::typeName)
    {
        const plasticityModel& plasticity = 
            refCast<const plasticityModel>(rheology);

	mu = plasticity.newMu().boundaryField()[patch().index()];
	lambda = plasticity.newLambda().boundaryField()[patch().index()];
    }


    //---------------------------//
    //- required fields
    //---------------------------//
    vectorField n = patch().nf();

    //- gradient of the field
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>("grad(" + fieldName_ + ")");



    //---------------------------//
    //- Set value fraction to fix reference patch normal
    //---------------------------//
    this->valueFraction() = sqr(n);

    if(nonLinear_ != OFF)
      {
        tensorField F = I + gradField;
        tensorField Finv = inv(F);
	scalarField J = det(F);
        vectorField nCurrent = Finv & n;
        nCurrent /= mag(nCurrent);
	this->valueFraction() = sqr(nCurrent);
      }

    //---------------------------//
    //- calculate the traction to apply
    //- set the shear to zero and leave the normal component
    //---------------------------//
    vectorField Traction(n.size(),vector::zero);

    //- incremental solvers
    if(fieldName_ == "DU")
      {
	const fvPatchField<symmTensor>& sigma =
	  patch().lookupPatchField<volSymmTensorField, symmTensor>("sigma");
	
	//- increment of traction
	Traction = -(n & sigma);
      }

    //---------------------------//
    //- calculate the normal gradient based on the traction
    //---------------------------//
    vectorField newGradient = 
      Traction
      - (n & (mu*gradField.T() - (mu + lambda)*gradField))
      - n*lambda*tr(gradField);

    //- if there is plasticity
    if(rheology.type() == plasticityModel::typeName)
    {
        const plasticityModel& plasticity = 
            refCast<const plasticityModel>(rheology);

        newGradient +=
	  2*mu*(n & plasticity.DEpsilonP().boundaryField()[patch().index()]);
    }

    //- if there are thermal effects
    if(this->db().objectRegistry::foundObject<thermalModel>("thermalProperties"))
      {
	const thermalModel& thermo =
	  this->db().objectRegistry::lookupObject<thermalModel>("thermalProperties");
    
	const fvPatchField<scalar>& T =
	  patch().lookupPatchField<volScalarField, scalar>("T");
      
	const scalarField threeKalpha =
	  (3*lambda + 2*mu)*
	  thermo.alpha()().boundaryField()[patch().index()];
      
	const scalarField T0 = thermo.T0()().boundaryField()[patch().index()];
	
	newGradient +=  (n*threeKalpha*(T - T0));
    }

    //- higher order non-linear terms
    if(nonLinear_ == UPDATED_LAGRANGIAN || nonLinear_ == TOTAL_LAGRANGIAN)
      {
	newGradient -=
	  (n & (mu*(gradField & gradField.T())))
	  + 0.5*n*lambda*(gradField && gradField);
	//- tensorial identity
	//- tr(gradField & gradField.T())*I == (gradField && gradField)*I
      }

    newGradient /= (2.0*mu + lambda);

    refGrad() = newGradient;

    directionMixedFvPatchVectorField::updateCoeffs();
}


// Write
void fixedDisplacementZeroShearFvPatchVectorField::write(Ostream& os) const
{
    directionMixedFvPatchVectorField::write(os);
    os.writeKeyword("nonLinear") << nonLinearNames_[nonLinear_] << token::END_STATEMENT << nl;
}


template<>
const char* Foam::NamedEnum<Foam::fixedDisplacementZeroShearFvPatchVectorField::nonLinearType, 3>::names[] =
  {
    "off",
    "updatedLagrangian",
    "totalLagrangian"
  };

const Foam::NamedEnum<Foam::fixedDisplacementZeroShearFvPatchVectorField::nonLinearType, 3>
Foam::fixedDisplacementZeroShearFvPatchVectorField::nonLinearNames_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, fixedDisplacementZeroShearFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
