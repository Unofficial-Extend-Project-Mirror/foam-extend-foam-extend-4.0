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

#include "fixedRotationFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
//#include "mathematicalConstants.H"
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
    rotationOrigin_(vector::zero),
    origCf_(0, vector::zero),
    nonLinear_(nonLinearGeometry::OFF)
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
    rotationOrigin_(ptf.rotationOrigin_),
    origCf_(ptf.origCf_),
    nonLinear_(nonLinearGeometry::OFF)
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
    rotationOrigin_(dict.lookup("rotationOrigin")),
    origCf_(patch().patch().faceCentres()),
    nonLinear_(nonLinearGeometry::OFF)
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
            Info<< "\tnonLinear set to updated Lagrangian"
                << endl;
        }
        else if (nonLinear_ == nonLinearGeometry::TOTAL_LAGRANGIAN)
        {
            Info << "\tnonLinear set to total Lagrangian"
                << endl;
        }
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


fixedRotationFvPatchVectorField::fixedRotationFvPatchVectorField
(
    const fixedRotationFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    fieldName_(pivpvf.fieldName_),
    rotationAngle_(pivpvf.rotationAngle_),
    rotationAxis_(pivpvf.rotationAxis_),
    rotationOrigin_(pivpvf.rotationOrigin_),
    origCf_(pivpvf.origCf_),
    nonLinear_(pivpvf.nonLinear_)
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
    rotationOrigin_(pivpvf.rotationOrigin_),
    origCf_(pivpvf.origCf_),
    nonLinear_(pivpvf.nonLinear_)
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

    tensor rotMat = RodriguesRotation(rotationAxis_, rotationAngle_);

    vectorField oldFaceCentres =
        dimensionedInternalField().mesh().C().boundaryField()[patch().index()];

    vectorField newFaceCentres =
        (rotMat & (oldFaceCentres - rotationOrigin_)) + rotationOrigin_;

    vectorField disp = newFaceCentres - oldFaceCentres;

    if (fieldName_ == "DU")
      {
        const fvPatchField<vector>& U =
          patch().lookupPatchField<volVectorField, vector>("U");
        disp -= U;
      }
    else if (fieldName_ != "U")
      {
        FatalError << "The displacement field should be U or DU"
                   << exit(FatalError);
      }

    fvPatchField<vector>::operator==
    (
        disp
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
    os.writeKeyword("rotationAngle")
        << rotationAngle_
        << token::END_STATEMENT << nl;
    os.writeKeyword("rotationAxis")
        << rotationAxis_
        << token::END_STATEMENT << nl;
    os.writeKeyword("rotationOrigin")
        << rotationOrigin_
        << token::END_STATEMENT << nl;
    os.writeKeyword("nonLinear")
        << nonLinearGeometry::nonLinearNames_[nonLinear_]
        << token::END_STATEMENT << nl;
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
