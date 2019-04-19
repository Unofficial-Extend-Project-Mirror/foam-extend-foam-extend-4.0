/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "swirlInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swirlInletVelocityFvPatchVectorField::
swirlInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    origin_(vector::zero),
    axis_(vector::zero),
    axialVelocity_(),
    radialVelocity_(),
    tangentialVelocity_()
{}


Foam::swirlInletVelocityFvPatchVectorField::
swirlInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),
    axialVelocity_(Function1<scalar>::New("axialVelocity", dict)),
    radialVelocity_(Function1<scalar>::New("radialVelocity", dict)),
    tangentialVelocity_(Function1<scalar>::New("tangentialVelocity", dict))
{}


Foam::swirlInletVelocityFvPatchVectorField::
swirlInletVelocityFvPatchVectorField
(
    const swirlInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    axialVelocity_(ptf.axialVelocity_, false),
    radialVelocity_(ptf.radialVelocity_, false),
    tangentialVelocity_(ptf.tangentialVelocity_, false)
{}


Foam::swirlInletVelocityFvPatchVectorField::
swirlInletVelocityFvPatchVectorField
(
    const swirlInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    axialVelocity_(ptf.axialVelocity_, false),
    radialVelocity_(ptf.radialVelocity_, false),
    tangentialVelocity_(ptf.tangentialVelocity_, false)
{}


Foam::swirlInletVelocityFvPatchVectorField::
swirlInletVelocityFvPatchVectorField
(
    const swirlInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    axialVelocity_(ptf.axialVelocity_, false),
    radialVelocity_(ptf.radialVelocity_, false),
    tangentialVelocity_(ptf.tangentialVelocity_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::swirlInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const vector axisHat = axis_/mag(axis_);

    // Radius vector in plane of rotation
    vectorField r(patch().Cf() - origin_);
    r -= (axisHat & r)*axisHat;
    const scalarField magr(mag(r));
    const vectorField rHat(r/magr);

    const scalarField axialVelocity(axialVelocity_->value(magr));
    const scalarField radialVelocity(radialVelocity_->value(magr));
    const scalarField tangentialVelocity(tangentialVelocity_->value(magr));

    operator==
    (
        axialVelocity*axisHat
      + radialVelocity*rHat
      + tangentialVelocity*(axisHat ^ rHat)
    );

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::swirlInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("origin") << origin_ << token::END_STATEMENT << nl;
    os.writeKeyword("axis") << axis_ << token::END_STATEMENT << nl;
    axialVelocity_->writeData(os);
    radialVelocity_->writeData(os);
    tangentialVelocity_->writeData(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       swirlInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
