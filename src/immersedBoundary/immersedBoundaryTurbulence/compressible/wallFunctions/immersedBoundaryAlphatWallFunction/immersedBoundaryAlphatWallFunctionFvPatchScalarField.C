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

#include "immersedBoundaryAlphatWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "standAlonePatch.H"
#include "surfaceWriter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

immersedBoundaryAlphatWallFunctionFvPatchScalarField::
immersedBoundaryAlphatWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatWallFunctionFvPatchScalarField(p, iF),
    immersedBoundaryFieldBase<scalar>(p, true, 1e-6)
{}


immersedBoundaryAlphatWallFunctionFvPatchScalarField::
immersedBoundaryAlphatWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphatWallFunctionFvPatchScalarField(p, iF),   // Do not read mixed data
    immersedBoundaryFieldBase<scalar>(p, true, 1e-6)
{
    this->readPatchType(dict);

    if (!isType<immersedBoundaryFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    scalarField::operator=(this->patchInternalField());
}


immersedBoundaryAlphatWallFunctionFvPatchScalarField::
immersedBoundaryAlphatWallFunctionFvPatchScalarField
(
    const immersedBoundaryAlphatWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphatWallFunctionFvPatchScalarField(p, iF), // Do not map mixed data.  Set patchType later
    immersedBoundaryFieldBase<scalar>(p, true, 1e-6)
{
    // Note: NO MAPPING.  Fields are created on the immersed boundary
    // HJ, 12/Apr/2012
    if (!isType<immersedBoundaryFvPatch>(p))
    {
        FatalErrorInFunction
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    this->setPatchType(ptf);

    // On creation of the mapped field, the internal field is dummy and
    // cannot be used.  Initialise the value to avoid errors
    // HJ, 1/Dec/2017
    scalarField::operator=(scalar(0));
}


immersedBoundaryAlphatWallFunctionFvPatchScalarField::
immersedBoundaryAlphatWallFunctionFvPatchScalarField
(
    const immersedBoundaryAlphatWallFunctionFvPatchScalarField& ptf
)
:
    alphatWallFunctionFvPatchScalarField(ptf),
    immersedBoundaryFieldBase<scalar>(ptf.ibPatch(), true, 1e-6)
{
    this->setPatchType(ptf);
}


immersedBoundaryAlphatWallFunctionFvPatchScalarField::
immersedBoundaryAlphatWallFunctionFvPatchScalarField
(
    const immersedBoundaryAlphatWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatWallFunctionFvPatchScalarField(ptf, iF),
    immersedBoundaryFieldBase<scalar>(ptf.ibPatch(), true, 1e-6)
{
    this->setPatchType(ptf);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void immersedBoundaryAlphatWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper&
)
{
    scalarField::operator=(this->patchInternalField());
}


void immersedBoundaryAlphatWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList&
)
{}


void immersedBoundaryAlphatWallFunctionFvPatchScalarField::updateOnMotion()
{
    if (size() != ibPatch().size())
    {
        // Use internal values, resizing the file if needed
        scalarField::operator=(this->patchInternalField());
    }
}


void immersedBoundaryAlphatWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Resize fields
    if (size() != patch().size())
    {
        Info<< "Resizing immersedBoundaryAlphatWallFunction in evaluate"
            << endl;

        scalarField::operator=(patchInternalField());
    }

    // Get non-constant reference to internal field
    scalarField& intField = const_cast<scalarField&>(this->internalField());

    // Set dead value
    this->setDeadValues(intField);

    alphatWallFunctionFvPatchScalarField::evaluate(commsType);
}


void immersedBoundaryAlphatWallFunctionFvPatchScalarField::write
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
    immersedBoundaryAlphatWallFunctionFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
