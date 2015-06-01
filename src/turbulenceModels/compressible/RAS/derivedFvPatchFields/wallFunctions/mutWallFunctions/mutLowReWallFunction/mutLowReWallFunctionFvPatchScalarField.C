/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "mutLowReWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> mutLowReWallFunctionFvPatchScalarField::calcMut() const
{
    return tmp<scalarField>(new scalarField(patch().size(), 0.0));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mutLowReWallFunctionFvPatchScalarField::mutLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutWallFunctionFvPatchScalarField(p, iF)
{}


mutLowReWallFunctionFvPatchScalarField::mutLowReWallFunctionFvPatchScalarField
(
    const mutLowReWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


mutLowReWallFunctionFvPatchScalarField::mutLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mutWallFunctionFvPatchScalarField(p, iF, dict)
{}


mutLowReWallFunctionFvPatchScalarField::mutLowReWallFunctionFvPatchScalarField
(
    const mutLowReWallFunctionFvPatchScalarField& mlrwfpsf
)
:
    mutWallFunctionFvPatchScalarField(mlrwfpsf)
{}


mutLowReWallFunctionFvPatchScalarField::mutLowReWallFunctionFvPatchScalarField
(
    const mutLowReWallFunctionFvPatchScalarField& mlrwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutWallFunctionFvPatchScalarField(mlrwfpsf, iF)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    mutLowReWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
