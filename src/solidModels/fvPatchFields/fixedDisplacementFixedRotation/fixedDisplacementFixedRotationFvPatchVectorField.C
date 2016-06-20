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

#include "fixedDisplacementFixedRotationFvPatchVectorField.H"
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

fixedDisplacementFixedRotationFvPatchVectorField::
fixedDisplacementFixedRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    fieldName_("undefined"),
    displacementTimeSeries_(),
    rotationAngleTimeSeries_(),
    rotationAxis_(vector::zero),
    rotationOrigin_(vector::zero)
{}


fixedDisplacementFixedRotationFvPatchVectorField::
fixedDisplacementFixedRotationFvPatchVectorField
(
    const fixedDisplacementFixedRotationFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_),
    displacementTimeSeries_(ptf.displacementTimeSeries_),
    rotationAngleTimeSeries_(ptf.rotationAngleTimeSeries_),
    rotationAxis_(ptf.rotationAxis_),
    rotationOrigin_(ptf.rotationOrigin_)
{}


fixedDisplacementFixedRotationFvPatchVectorField::
fixedDisplacementFixedRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    fieldName_(dimensionedInternalField().name()),
    displacementTimeSeries_(dict.subDict("timeVaryingDisplacement")),
    rotationAngleTimeSeries_(dict.subDict("timeVaryingRotation")),
    rotationAxis_(dict.subDict("timeVaryingRotation").lookup("rotationAxis")),
    rotationOrigin_
    (dict.subDict("timeVaryingRotation").lookup("rotationOrigin"))
{
  Info << "fixedDisplacementFixedRotation boundary condition"
       << endl;

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


fixedDisplacementFixedRotationFvPatchVectorField::
fixedDisplacementFixedRotationFvPatchVectorField
(
    const fixedDisplacementFixedRotationFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    fieldName_(pivpvf.fieldName_),
    displacementTimeSeries_(pivpvf.displacementTimeSeries_),
    rotationAngleTimeSeries_(pivpvf.rotationAngleTimeSeries_),
    rotationAxis_(pivpvf.rotationAxis_),
    rotationOrigin_(pivpvf.rotationOrigin_)
{}


fixedDisplacementFixedRotationFvPatchVectorField::
fixedDisplacementFixedRotationFvPatchVectorField
(
    const fixedDisplacementFixedRotationFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    fieldName_(pivpvf.fieldName_),
    displacementTimeSeries_(pivpvf.displacementTimeSeries_),
    rotationAngleTimeSeries_(pivpvf.rotationAngleTimeSeries_),
    rotationAxis_(pivpvf.rotationAxis_),
    rotationOrigin_(pivpvf.rotationOrigin_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedDisplacementFixedRotationFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Patch displacement is a superposition of a translation and a rotation


    // Rotation

    // Look up angle depending on current time
    scalar angle =
        rotationAngleTimeSeries_(this->db().time().timeOutputValue());

    tensor rotMat = RodriguesRotation(rotationAxis_, angle);

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
        FatalError
            << "The displacement field should be U or DU"
            << exit(FatalError);
      }

    // translation
    disp += displacementTimeSeries_(this->db().time().timeOutputValue());

    fvPatchField<vector>::operator==
      (
       disp
       );

    fixedValueFvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::Field<vector> >
fixedDisplacementFixedRotationFvPatchVectorField::
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

tmp<Field<vector> > fixedDisplacementFixedRotationFvPatchVectorField::
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

void fixedDisplacementFixedRotationFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);

    os.writeKeyword("timeVaryingDisplacement") << nl;
    os << token::BEGIN_BLOCK << nl;
    displacementTimeSeries_.write(os);
    os << token::END_BLOCK << nl;

    os.writeKeyword("timeVaryingRotation") << nl;
    os << token::BEGIN_BLOCK << nl;
    rotationAngleTimeSeries_.write(os);
    os.writeKeyword("rotationAxis")
        << rotationAxis_ << token::END_STATEMENT << nl;
    os.writeKeyword("rotationOrigin")
        << rotationOrigin_ << token::END_STATEMENT << nl;
    os << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedDisplacementFixedRotationFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
