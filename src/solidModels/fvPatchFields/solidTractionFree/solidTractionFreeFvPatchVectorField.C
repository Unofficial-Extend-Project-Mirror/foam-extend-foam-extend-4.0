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
//#include "rheologyModel.H"
//#include "plasticityModel.H"
//#include "thermalModel.H"
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
    nonLinear_(OFF),
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
    nonLinear_(OFF),
    orthotropic_(false)
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

    gradient() = tractionBoundaryGradient()
      (
       vectorField(patch().size(), vector::zero),
       scalarField(patch().size(), 0.0),
       word(fieldName_),
       patch(),
       orthotropic_,
       NamedEnum<Foam::solidTractionFreeFvPatchVectorField::nonLinearType, 3>::names[nonLinear_]
       )();

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
