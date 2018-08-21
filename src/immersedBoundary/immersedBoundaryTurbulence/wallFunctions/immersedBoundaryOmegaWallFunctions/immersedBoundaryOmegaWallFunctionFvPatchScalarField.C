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

#include "immersedBoundaryOmegaWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

immersedBoundaryOmegaWallFunctionFvPatchScalarField::
immersedBoundaryOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    omegaWallFunctionFvPatchScalarField(p, iF),
    immersedBoundaryFieldBase<scalar>(p, true, 90)
{}


immersedBoundaryOmegaWallFunctionFvPatchScalarField::
immersedBoundaryOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    omegaWallFunctionFvPatchScalarField(p, iF, dict),
    immersedBoundaryFieldBase<scalar>(p, true, 90)
{}


immersedBoundaryOmegaWallFunctionFvPatchScalarField::
immersedBoundaryOmegaWallFunctionFvPatchScalarField
(
    const immersedBoundaryOmegaWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    omegaWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    immersedBoundaryFieldBase<scalar>(p, true, 90)
{}


immersedBoundaryOmegaWallFunctionFvPatchScalarField::
immersedBoundaryOmegaWallFunctionFvPatchScalarField
(
    const immersedBoundaryOmegaWallFunctionFvPatchScalarField& ptf
)
:
    omegaWallFunctionFvPatchScalarField(ptf),
    immersedBoundaryFieldBase<scalar>(ptf.ibPatch(), true, 90)
{}


immersedBoundaryOmegaWallFunctionFvPatchScalarField::
immersedBoundaryOmegaWallFunctionFvPatchScalarField
(
    const immersedBoundaryOmegaWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    omegaWallFunctionFvPatchScalarField(ptf, iF),
    immersedBoundaryFieldBase<scalar>(ptf.ibPatch(), true, 90)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void immersedBoundaryOmegaWallFunctionFvPatchScalarField::updateOnMotion()
{
    if (size() != ibPatch().size())
    {
        // Use internal values, resizing the file if needed
        scalarField::operator=(this->patchInternalField());

        // Resize refValue as well.  HJ, 10/Jul/2018
        refValue() = this->patchInternalField();
    }
}


void immersedBoundaryOmegaWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Resize fields
    if (size() != patch().size())
    {
        Info<< "Resizing immersedBoundaryEpsilonWallFunction in evaluate: "
            << "patch: " << patch().size() << " field: " << size()
            << endl;

        *this == patchInternalField();
        refValue() = patchInternalField();
    }

    // If G field is present, execute evaluation
    // Remove the warning from the IB patch
    // HJ, 20/May/2018
    if (db().foundObject<volScalarField>(GName()))
    {
        omegaWallFunctionFvPatchScalarField::updateCoeffs();
    }
}


void immersedBoundaryOmegaWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeLocalEntries(os);

    // The value entry needs to be written with zero size
    scalarField::null().writeEntry("value", os);
    // this->writeEntry("value", os);

    writeField(*this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    immersedBoundaryOmegaWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
