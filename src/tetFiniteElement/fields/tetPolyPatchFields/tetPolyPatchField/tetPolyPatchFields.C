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

Description

\*---------------------------------------------------------------------------*/

#include "tetFemMatrices.H"
#include "tetPolyPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

typedef PointPatchField
<tetPolyPatchField, tetPointMesh, tetPolyPatch, tetFemMatrix, scalar>
    baseTetPolyPatchScalarField;

typedef PointPatchField
<tetPolyPatchField, tetPointMesh, tetPolyPatch, tetFemMatrix, vector>
    baseTetPolyPatchVectorField;

typedef PointPatchField
<tetPolyPatchField, tetPointMesh, tetPolyPatch, tetFemMatrix, sphericalTensor>
    baseTetPolyPatchSphericalTensorField;

typedef PointPatchField
<tetPolyPatchField, tetPointMesh, tetPolyPatch, tetFemMatrix, symmTensor>
    baseTetPolyPatchSymmTensorField;

typedef PointPatchField
<tetPolyPatchField, tetPointMesh, tetPolyPatch, tetFemMatrix, symmTensor4thOrder>
    baseTetPolyPatchSymmTensor4thOrderField;

typedef PointPatchField
<tetPolyPatchField, tetPointMesh, tetPolyPatch, tetFemMatrix, diagTensor>
    baseTetPolyPatchDiagTensorField;

typedef PointPatchField
<tetPolyPatchField, tetPointMesh, tetPolyPatch, tetFemMatrix, tensor>
    baseTetPolyPatchTensorField;

defineNamedTemplateTypeNameAndDebug(baseTetPolyPatchScalarField, 0);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchScalarField, PointPatch);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchScalarField, patchMapper);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchScalarField, dictionary);

defineNamedTemplateTypeNameAndDebug(baseTetPolyPatchVectorField, 0);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchVectorField, PointPatch);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchVectorField, patchMapper);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchVectorField, dictionary);

defineNamedTemplateTypeNameAndDebug(baseTetPolyPatchSphericalTensorField, 0);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchSphericalTensorField, PointPatch);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchSphericalTensorField, patchMapper);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchSphericalTensorField, dictionary);

defineNamedTemplateTypeNameAndDebug(baseTetPolyPatchSymmTensorField, 0);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchSymmTensorField, PointPatch);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchSymmTensorField, patchMapper);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchSymmTensorField, dictionary);

defineNamedTemplateTypeNameAndDebug(baseTetPolyPatchSymmTensor4thOrderField, 0);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchSymmTensor4thOrderField, PointPatch);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchSymmTensor4thOrderField, patchMapper);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchSymmTensor4thOrderField, dictionary);

defineNamedTemplateTypeNameAndDebug(baseTetPolyPatchDiagTensorField, 0);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchDiagTensorField, PointPatch);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchDiagTensorField, patchMapper);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchDiagTensorField, dictionary);

defineNamedTemplateTypeNameAndDebug(baseTetPolyPatchTensorField, 0);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchTensorField, PointPatch);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchTensorField, patchMapper);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchTensorField, dictionary);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
