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
    solidDirectionMixedFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "solidDirectionMixedFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidDirectionMixedFvPatchVectorField::solidDirectionMixedFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    fieldName_("undefined")
{}


solidDirectionMixedFvPatchVectorField::solidDirectionMixedFvPatchVectorField
(
    const solidDirectionMixedFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_)
{}


solidDirectionMixedFvPatchVectorField::solidDirectionMixedFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
  directionMixedFvPatchVectorField(p, iF, dict),
  fieldName_(dimensionedInternalField().name())
{
  Field<vector> normalValue = transform(valueFraction(), refValue());

  Field<vector> gradValue =
  this->patchInternalField() + refGrad()/this->patch().deltaCoeffs();

  //- non-ortho corrected gradValue
  //- gradField will not have been created so I must do this during updateCoeffs
  /*const fvPatchField<tensor>& gradField =
    patch().lookupPatchField<volTensorField, tensor>("grad(" +fieldName_ + ")");
    vectorField n = patch().nf();
    vectorField delta = patch().delta();
    vectorField k = delta - n*(n&delta);
    Field<vector> gradValue = this->patchInternalField()
    + (k&gradField.patchInternalField())
    + refGrad()/this->patch().deltaCoeffs();
  */

  Field<vector> transformGradValue =
    transform(I - valueFraction(), gradValue);

  Field<vector>::operator=(normalValue + transformGradValue);
}

solidDirectionMixedFvPatchVectorField::solidDirectionMixedFvPatchVectorField
(
    const solidDirectionMixedFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF),
    fieldName_(ptf.fieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void solidDirectionMixedFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void solidDirectionMixedFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);
}


void solidDirectionMixedFvPatchVectorField::updateCoeffs()
{
    directionMixedFvPatchVectorField::updateCoeffs();
}

void solidDirectionMixedFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<vector> normalValue = transform(valueFraction(), refValue());

    //- no correction
    //Field<vector> gradValue =
    //this->patchInternalField() + refGrad()/this->patch().deltaCoeffs();

    //- non-ortho corrected gradValue
    const fvPatchField<tensor>& gradField =
      patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" +fieldName_ + ")"
            );
    vectorField n = patch().nf();
    vectorField delta = patch().delta();
    vectorField k = delta - n*(n&delta);
    Field<vector> gradValue = this->patchInternalField()
      + (k&gradField.patchInternalField())
      + refGrad()/this->patch().deltaCoeffs();

    Field<vector> transformGradValue =
      transform(I - valueFraction(), gradValue);

    Field<vector>::operator=(normalValue + transformGradValue);

    fvPatchField<vector>::evaluate();
}

Foam::tmp<Foam::Field<vector> >
solidDirectionMixedFvPatchVectorField::snGrad() const
{
    Field<vector> pif = this->patchInternalField();

  //- fixedValue snGrad with no correction
  //  return (*this - patchInternalField())*this->patch().deltaCoeffs();
      Field<vector> normalValue = transform(valueFraction(), refValue());

    const fvPatchField<tensor>& gradField =
      patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" +fieldName_ + ")"
            );
    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();
    //- correction vector
    vectorField k = delta - n*(n&delta);

    Field<vector> gradValue =
      (this->refGrad()/this->patch().deltaCoeffs())
      + (pif)
      + (k&gradField.patchInternalField());

    Field<vector> transformGradValue =
      transform(I - (this->valueFraction()), gradValue);

    Field<vector> patchValue = normalValue + transformGradValue;

    return
      (
       patchValue
       - ((this->patchInternalField()) + (k&gradField.patchInternalField()))
       )*this->patch().deltaCoeffs();
}

// Write
void solidDirectionMixedFvPatchVectorField::write(Ostream& os) const
{
    directionMixedFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidDirectionMixedFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
