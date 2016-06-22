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

#include "pressureInletUniformVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pressureInletUniformVelocityFvPatchVectorField::
pressureInletUniformVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    pressureInletVelocityFvPatchVectorField(p, iF)
{}


pressureInletUniformVelocityFvPatchVectorField::
pressureInletUniformVelocityFvPatchVectorField
(
    const pressureInletUniformVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    pressureInletVelocityFvPatchVectorField(ptf, p, iF, mapper)
{}


pressureInletUniformVelocityFvPatchVectorField::
pressureInletUniformVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    pressureInletVelocityFvPatchVectorField(p, iF, dict)
{}


pressureInletUniformVelocityFvPatchVectorField::
pressureInletUniformVelocityFvPatchVectorField
(
    const pressureInletUniformVelocityFvPatchVectorField& pivpvf
)
:
    pressureInletVelocityFvPatchVectorField(pivpvf)
{}


pressureInletUniformVelocityFvPatchVectorField::
pressureInletUniformVelocityFvPatchVectorField
(
    const pressureInletUniformVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    pressureInletVelocityFvPatchVectorField(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pressureInletUniformVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    pressureInletVelocityFvPatchVectorField::updateCoeffs();

    operator==(patch().nf()*gSum(patch().Sf() & *this)/gSum(patch().magSf()));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void pressureInletUniformVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    operator==(patch().nf()*gSum(patch().Sf() & pvf)/gSum(patch().magSf()));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    pressureInletUniformVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
