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

#include "nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField.H"
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

nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField::
nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutURoughWallFunctionFvPatchScalarField(p, iF)
{}


nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField::
nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
(
    const nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutURoughWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField::
nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutURoughWallFunctionFvPatchScalarField(p, iF, dict)
{}


nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField::
nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
(
    const nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField& rwfpsf
)
:
    nutURoughWallFunctionFvPatchScalarField(rwfpsf)
{}


nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField::
nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
(
    const nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutURoughWallFunctionFvPatchScalarField(rwfpsf, iF)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
