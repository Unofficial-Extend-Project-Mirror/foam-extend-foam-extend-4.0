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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fixedRotationFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "mathematicalConstants.H"
#include "RodriguesRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedRotationFvPatchVectorField::fixedRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    fieldName_("undefined"),
    rotationAngle_(0.0),
    rotationAxis_(vector::zero),
    rotationOrigin_(vector::zero)
{}


fixedRotationFvPatchVectorField::fixedRotationFvPatchVectorField
(
    const fixedRotationFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_),
    rotationAngle_(ptf.rotationAngle_),
    rotationAxis_(ptf.rotationAxis_),
    rotationOrigin_(ptf.rotationOrigin_)
{}


fixedRotationFvPatchVectorField::fixedRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    fieldName_(dimensionedInternalField().name()),
    rotationAngle_(readScalar(dict.lookup("rotationAngle"))),
    rotationAxis_(dict.lookup("rotationAxis")),
    rotationOrigin_(dict.lookup("rotationOrigin"))
{
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


fixedRotationFvPatchVectorField::fixedRotationFvPatchVectorField
(
    const fixedRotationFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    fieldName_(pivpvf.fieldName_),
    rotationAngle_(pivpvf.rotationAngle_),
    rotationAxis_(pivpvf.rotationAxis_),
    rotationOrigin_(pivpvf.rotationOrigin_)
{}


fixedRotationFvPatchVectorField::fixedRotationFvPatchVectorField
(
    const fixedRotationFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    fieldName_(pivpvf.fieldName_),
    rotationAngle_(pivpvf.rotationAngle_),
    rotationAxis_(pivpvf.rotationAxis_),
    rotationOrigin_(pivpvf.rotationOrigin_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<vector> > fixedRotationFvPatchVectorField::
snGrad() const
{
  //- fixedValue snGrad with no correction
  //  return (*this - patchInternalField())*this->patch().deltaCoeffs();
  
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" +fieldName_ + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();

    //- correction vector
    vectorField k = delta - n*(n&delta);

    return 
    (
        *this 
      - (patchInternalField() + (k&gradField.patchInternalField()))
      )*this->patch().deltaCoeffs();
}

void fixedRotationFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //- convert rotation to radians
    //    scalar theta = rotationAngle_ *  Foam::mathematicalConstant::pi / 180;

    //- create rotation matrix
    // tensor rotMat(::cos(theta), -(::sin(theta)), 0,
    // 		  ::sin(theta), ::cos(theta),    0,
    // 		  0,            0,               1);

    tensor rotMat = RodriguesRotation(rotationAxis_, rotationAngle_);

    const vectorField& oldFaceCentres = dimensionedInternalField().mesh().C().boundaryField()[patch().index()];
	  
    vectorField newFaceCentres = (rotMat & (oldFaceCentres - rotationOrigin_)) + rotationOrigin_;
    
    fvPatchField<vector>::operator==
    (
        newFaceCentres - oldFaceCentres
    );
    fixedValueFvPatchVectorField::updateCoeffs();
}


tmp<Field<vector> > fixedRotationFvPatchVectorField::
gradientBoundaryCoeffs() const
{
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" +fieldName_ + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();

    //- correction vector
    vectorField k = delta - n*(n&delta);

    return this->patch().deltaCoeffs()
       *(*this - (k&gradField.patchInternalField()));
}

void fixedRotationFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);
    os.writeKeyword("rotationAngle") << rotationAngle_ << token::END_STATEMENT << nl;
    os.writeKeyword("rotationAxis") << rotationAxis_ << token::END_STATEMENT << nl;
    os.writeKeyword("rotationOrigin") << rotationOrigin_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedRotationFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
