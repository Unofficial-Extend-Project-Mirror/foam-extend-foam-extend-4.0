/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

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

defineNamedTemplateTypeNameAndDebug(baseTetPolyPatchTensorField, 0);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchTensorField, PointPatch);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchTensorField, patchMapper);
defineTemplateRunTimeSelectionTable(baseTetPolyPatchTensorField, dictionary);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
