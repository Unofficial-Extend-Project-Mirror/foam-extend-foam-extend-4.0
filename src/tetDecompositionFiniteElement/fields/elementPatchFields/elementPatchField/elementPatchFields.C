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
#include "elementPatchFields.H"
#include "elementMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

typedef PointPatchField
<elementPatchField, elementMesh, tetPolyPatch, DummyMatrix, scalar>
    baseElementPatchScalarField;

typedef PointPatchField
<elementPatchField, elementMesh, tetPolyPatch, DummyMatrix, vector>
    baseElementPatchVectorField;

typedef PointPatchField
<elementPatchField, elementMesh, tetPolyPatch, DummyMatrix, sphericalTensor>
    baseElementPatchSphericalTensorField;

typedef PointPatchField
<elementPatchField, elementMesh, tetPolyPatch, DummyMatrix, symmTensor>
    baseElementPatchSymmTensorField;

typedef PointPatchField
<elementPatchField, elementMesh, tetPolyPatch, DummyMatrix, tensor>
    baseElementPatchTensorField;

defineNamedTemplateTypeNameAndDebug(baseElementPatchScalarField, 0);
defineTemplateRunTimeSelectionTable(baseElementPatchScalarField, PointPatch);
defineTemplateRunTimeSelectionTable(baseElementPatchScalarField, patchMapper);
defineTemplateRunTimeSelectionTable(baseElementPatchScalarField, dictionary);

defineNamedTemplateTypeNameAndDebug(baseElementPatchVectorField, 0);
defineTemplateRunTimeSelectionTable(baseElementPatchVectorField, PointPatch);
defineTemplateRunTimeSelectionTable(baseElementPatchVectorField, patchMapper);
defineTemplateRunTimeSelectionTable(baseElementPatchVectorField, dictionary);

defineNamedTemplateTypeNameAndDebug(baseElementPatchSphericalTensorField, 0);
defineTemplateRunTimeSelectionTable(baseElementPatchSphericalTensorField, PointPatch);
defineTemplateRunTimeSelectionTable(baseElementPatchSphericalTensorField, patchMapper);
defineTemplateRunTimeSelectionTable(baseElementPatchSphericalTensorField, dictionary);

defineNamedTemplateTypeNameAndDebug(baseElementPatchSymmTensorField, 0);
defineTemplateRunTimeSelectionTable(baseElementPatchSymmTensorField, PointPatch);
defineTemplateRunTimeSelectionTable(baseElementPatchSymmTensorField, patchMapper);
defineTemplateRunTimeSelectionTable(baseElementPatchSymmTensorField, dictionary);

defineNamedTemplateTypeNameAndDebug(baseElementPatchTensorField, 0);
defineTemplateRunTimeSelectionTable(baseElementPatchTensorField, PointPatch);
defineTemplateRunTimeSelectionTable(baseElementPatchTensorField, patchMapper);
defineTemplateRunTimeSelectionTable(baseElementPatchTensorField, dictionary);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
