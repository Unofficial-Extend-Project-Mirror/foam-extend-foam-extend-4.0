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

#include "immersedBoundaryNutWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "standAlonePatch.H"
#include "surfaceWriter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

immersedBoundaryNutWallFunctionFvPatchScalarField::
immersedBoundaryNutWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(p, iF),
    immersedBoundaryFieldBase<scalar>(p, true, 1e-6)
{}


immersedBoundaryNutWallFunctionFvPatchScalarField::
immersedBoundaryNutWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutkWallFunctionFvPatchScalarField(p, iF, dict),
    immersedBoundaryFieldBase<scalar>(p, true, 1e-6)
{
    this->readPatchType(dict);
}


immersedBoundaryNutWallFunctionFvPatchScalarField::
immersedBoundaryNutWallFunctionFvPatchScalarField
(
    const immersedBoundaryNutWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    immersedBoundaryFieldBase<scalar>(p, true, 1e-6)
{
    this->setPatchType(ptf);
}


immersedBoundaryNutWallFunctionFvPatchScalarField::
immersedBoundaryNutWallFunctionFvPatchScalarField
(
    const immersedBoundaryNutWallFunctionFvPatchScalarField& ptf
)
:
    nutkWallFunctionFvPatchScalarField(ptf),
    immersedBoundaryFieldBase<scalar>(ptf.ibPatch(), true, 1e-6)
{
    this->setPatchType(ptf);
}


immersedBoundaryNutWallFunctionFvPatchScalarField::
immersedBoundaryNutWallFunctionFvPatchScalarField
(
    const immersedBoundaryNutWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(ptf, iF),
    immersedBoundaryFieldBase<scalar>(ptf.ibPatch(), true, 1e-6)
{
    this->setPatchType(ptf);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void immersedBoundaryNutWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper&
)
{
    scalarField::operator=(this->patchInternalField());
}


void immersedBoundaryNutWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList&
)
{}


void immersedBoundaryNutWallFunctionFvPatchScalarField::updateOnMotion()
{
    if (size() != ibPatch().size())
    {
        // Use internal values, resizing the file if needed
        scalarField::operator=(this->patchInternalField());
    }
}


void immersedBoundaryNutWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Resize fields
    if (size() != patch().size())
    {
        Info<< "Resizing immersedBoundaryNutWallFunction in evaluate"
            << endl;

        *this == patchInternalField();
    }

    // Get non-constant reference to internal field
    scalarField& intField = const_cast<scalarField&>(this->internalField());

    // Set dead value
    this->setDeadValues(intField);

    nutkWallFunctionFvPatchScalarField::evaluate(commsType);
}


void immersedBoundaryNutWallFunctionFvPatchScalarField::write(Ostream& os) const
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
    immersedBoundaryNutWallFunctionFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
