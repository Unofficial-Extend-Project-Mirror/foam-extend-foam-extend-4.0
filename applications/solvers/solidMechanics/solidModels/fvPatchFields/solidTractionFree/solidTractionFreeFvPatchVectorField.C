/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

\*---------------------------------------------------------------------------*/

#include "solidTractionFreeFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "rheologyModel.H"
#include "plasticityModel.H"
#include "thermalModel.H"

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
    nonLinear_(OFF)
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
    nonLinear_(OFF)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;

    Info << "Patch " << patch().name()
	 << "\tTraction boundary field: " << fieldName_ << endl;

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
    nonLinear_(stpvf.nonLinear_)
{}


solidTractionFreeFvPatchVectorField::
solidTractionFreeFvPatchVectorField
(
    const solidTractionFreeFvPatchVectorField& stpvf
)
:
    fixedGradientFvPatchVectorField(stpvf),
    fieldName_(stpvf.fieldName_),
    nonLinear_(stpvf.nonLinear_)
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
    nonLinear_(stpvf.nonLinear_)
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
    //- calculate the traction to apply
    //---------------------------//
    vectorField Traction(n.size(),vector::zero);

    //- incremental solvers
    if(fieldName_ == "DU")
      {
	const fvPatchField<symmTensor>& sigma =
	  patch().lookupPatchField<volSymmTensorField, symmTensor>("sigma");
	
	//- increment of traction
	Traction = - (n & sigma);
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
      
        const fvPatchField<scalar>& threeKalpha =
          patch().lookupPatchField<volScalarField, scalar>("((threeK*rho)*alpha)");
      
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

    gradient() = newGradient;

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
      + (k&gradField.patchInternalField())
      + gradient()/this->patch().deltaCoeffs()
    );

    fvPatchField<vector>::evaluate();
}

// Write
void solidTractionFreeFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("nonLinear") << nonLinearNames_[nonLinear_] << token::END_STATEMENT << nl;
    writeEntry("value", os);
}

template<>
const char* Foam::NamedEnum<Foam::solidTractionFreeFvPatchVectorField::nonLinearType, 3>::names[] =
  {
    "off",
    "updatedLagrangian",
    "totalLagrangian"
  };

const Foam::NamedEnum<Foam::solidTractionFreeFvPatchVectorField::nonLinearType, 3>
Foam::solidTractionFreeFvPatchVectorField::nonLinearNames_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidTractionFreeFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
