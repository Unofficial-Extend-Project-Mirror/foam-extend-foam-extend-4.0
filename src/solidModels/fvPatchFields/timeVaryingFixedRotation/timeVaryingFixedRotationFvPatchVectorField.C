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

#include "timeVaryingFixedRotationFvPatchVectorField.H"
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

timeVaryingFixedRotationFvPatchVectorField::
timeVaryingFixedRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    timeSeries_(),
    fieldName_("undefined"),
    origFaceCentres_(0.0,vector::zero),
    //rotationAngle_(0.0),
    rotationAxis_(vector::zero),
    rotationOrigin_(vector::zero)
{}



timeVaryingFixedRotationFvPatchVectorField::
timeVaryingFixedRotationFvPatchVectorField
(
    const timeVaryingFixedRotationFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    timeSeries_(ptf.timeSeries_),
    fieldName_(ptf.fieldName_),
    accumulatedAngle_(ptf.accumulatedAngle_),
    origFaceCentres_(ptf.origFaceCentres_),
    //rotationAngle_(ptf.rotationAngle_),
    rotationAxis_(ptf.rotationAxis_),
    rotationOrigin_(ptf.rotationOrigin_)
{}


timeVaryingFixedRotationFvPatchVectorField::
timeVaryingFixedRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    timeSeries_(dict),
    fieldName_(dimensionedInternalField().name()),
    accumulatedAngle_(0.0),
    origFaceCentres_(patch().patch().faceCentres()),
    //rotationAngle_(0.0), //readScalar(dict.lookup("rotationAngle"))),
    rotationAxis_(dict.lookup("rotationAxis")),
    rotationOrigin_(dict.lookup("rotationOrigin"))
{
  Info << "Patch " << patch().name() << " is timeVaryingFixedRotation" << endl;

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


timeVaryingFixedRotationFvPatchVectorField::
timeVaryingFixedRotationFvPatchVectorField
(
    const timeVaryingFixedRotationFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    timeSeries_(pivpvf.timeSeries_),
    fieldName_(pivpvf.fieldName_),
    accumulatedAngle_(pivpvf.accumulatedAngle_),
    origFaceCentres_(pivpvf.origFaceCentres_),
    //rotationAngle_(pivpvf.rotationAngle_),
    rotationAxis_(pivpvf.rotationAxis_),
    rotationOrigin_(pivpvf.rotationOrigin_)
{}


timeVaryingFixedRotationFvPatchVectorField::
timeVaryingFixedRotationFvPatchVectorField
(
    const timeVaryingFixedRotationFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    timeSeries_(pivpvf.timeSeries_),
    fieldName_(pivpvf.fieldName_),
    accumulatedAngle_(pivpvf.accumulatedAngle_),
    origFaceCentres_(pivpvf.origFaceCentres_),
    //rotationAngle_(pivpvf.rotationAngle_),
    rotationAxis_(pivpvf.rotationAxis_),
    rotationOrigin_(pivpvf.rotationOrigin_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<vector> > timeVaryingFixedRotationFvPatchVectorField::
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

void timeVaryingFixedRotationFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //- Look up angle depending on current time
    scalar angle = timeSeries_(this->db().time().timeOutputValue());

    //- time varying rotation
    tensor rotMat = RodriguesRotation(rotationAxis_, angle);

    vectorField newFaceCentres =
      (rotMat & (origFaceCentres_ - rotationOrigin_)) + rotationOrigin_;

    vectorField disp = newFaceCentres - origFaceCentres_;

    if (fieldName_ == "DU" )
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

    fvPatchField<vector>::operator==
    (
        disp
    );
    fixedValueFvPatchVectorField::updateCoeffs();
}


tmp<Field<vector> > timeVaryingFixedRotationFvPatchVectorField::
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

void timeVaryingFixedRotationFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);
    timeSeries_.write(os);
    os.writeKeyword("rotationAxis")
        << rotationAxis_ << token::END_STATEMENT << nl;
    os.writeKeyword("rotationOrigin")
        << rotationOrigin_ << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    timeVaryingFixedRotationFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
