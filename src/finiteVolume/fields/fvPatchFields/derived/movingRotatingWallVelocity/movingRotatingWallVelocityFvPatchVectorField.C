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

#include "movingRotatingWallVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingRotatingWallVelocityFvPatchVectorField::
movingRotatingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    centre_(0, 0, 0),
    axis_(1, 0, 0),
    rpm_(0)
{}


Foam::movingRotatingWallVelocityFvPatchVectorField::
movingRotatingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const point& centre,
    const vector& axis,
    const scalar rpm
)
:
    fixedValueFvPatchVectorField(p, iF),
    centre_(centre),
    axis_(axis),
    rpm_(rpm)
{
    if (mag(axis_) < SMALL)
    {
        FatalErrorIn
        (
            "movingRotatingWallVelocityFvPatchVectorField::"
            "movingRotatingWallVelocityFvPatchVectorField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<vector, volMesh>& iF,\n"
            "    const point& centre,\n"
            "    const vector& axis,\n"
            "    const scalar rpm\n"
            ")"
        )   << "Badly defined axis: zero magnitude: " << axis_
            << " for patch " << patch().name()
            << abort(FatalError);
    }

    axis_ /= mag(axis_);
}


Foam::movingRotatingWallVelocityFvPatchVectorField::
movingRotatingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    centre_(dict.lookup("centre")),
    axis_(dict.lookup("axis")),
    rpm_(readScalar(dict.lookup("rpm")))

{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    if (mag(axis_) < SMALL)
    {
        FatalErrorIn
        (
            "movingRotatingWallVelocityFvPatchVectorField::"
            "movingRotatingWallVelocityFvPatchVectorField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<vector, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Badly defined axis: zero magnitude: " << axis_
            << " for patch " << patch().name()
            << abort(FatalError);
    }

    axis_ /= mag(axis_);
}


Foam::movingRotatingWallVelocityFvPatchVectorField::
movingRotatingWallVelocityFvPatchVectorField
(
    const movingRotatingWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    centre_(ptf.centre_),
    axis_(ptf.axis_),
    rpm_(ptf.rpm_)
{}


Foam::movingRotatingWallVelocityFvPatchVectorField::
movingRotatingWallVelocityFvPatchVectorField
(
    const movingRotatingWallVelocityFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    centre_(pivpvf.centre_),
    axis_(pivpvf.axis_),
    rpm_(pivpvf.rpm_)
{}


Foam::movingRotatingWallVelocityFvPatchVectorField::
movingRotatingWallVelocityFvPatchVectorField
(
    const movingRotatingWallVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    centre_(pivpvf.centre_),
    axis_(pivpvf.axis_),
    rpm_(pivpvf.rpm_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingRotatingWallVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatch& p = patch();

    const volVectorField& U =
        db().lookupObject<volVectorField>(dimensionedInternalField().name());
    scalarField phip =
        p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U));

    vectorField n = p.nf();
    const scalarField& magSf = p.magSf();
    scalarField Un = phip/(magSf + VSMALL);

    vectorField Urot =
        (axis_ ^ (patch().Cf() - centre_))*rpm_*2*mathematicalConstant::pi/60;

    vectorField::operator=(Urot + n*(Un - (n & Urot)));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::movingRotatingWallVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("centre") << centre_ << token::END_STATEMENT << nl;
    os.writeKeyword("axis") << axis_ << token::END_STATEMENT << endl;
    os.writeKeyword("rpm") << rpm_ << token::END_STATEMENT << endl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        movingRotatingWallVelocityFvPatchVectorField
    );
} // End namespace Foam

// ************************************************************************* //
