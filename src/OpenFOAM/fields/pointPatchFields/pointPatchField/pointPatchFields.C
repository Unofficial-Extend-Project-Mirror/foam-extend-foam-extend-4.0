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

\*---------------------------------------------------------------------------*/

#include "pointPatchFields.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

typedef PointPatchField
<pointPatchField, pointMesh, pointPatch, DummyMatrix, scalar>
    basePointPatchScalarField;

typedef PointPatchField
<pointPatchField, pointMesh, pointPatch, DummyMatrix, vector>
    basePointPatchVectorField;

typedef PointPatchField
<pointPatchField, pointMesh, pointPatch, DummyMatrix, sphericalTensor>
    basePointPatchSphericalTensorField;

typedef PointPatchField
<pointPatchField, pointMesh, pointPatch, DummyMatrix, symmTensor>
    basePointPatchSymmTensorField;

typedef PointPatchField
<pointPatchField, pointMesh, pointPatch, DummyMatrix, tensor>
    basePointPatchTensorField;

defineNamedTemplateTypeNameAndDebug(basePointPatchScalarField, 0);
defineTemplateRunTimeSelectionTable(basePointPatchScalarField, PointPatch);
defineTemplateRunTimeSelectionTable(basePointPatchScalarField, patchMapper);
defineTemplateRunTimeSelectionTable(basePointPatchScalarField, dictionary);

defineNamedTemplateTypeNameAndDebug(basePointPatchVectorField, 0);
defineTemplateRunTimeSelectionTable(basePointPatchVectorField, PointPatch);
defineTemplateRunTimeSelectionTable(basePointPatchVectorField, patchMapper);
defineTemplateRunTimeSelectionTable(basePointPatchVectorField, dictionary);

defineNamedTemplateTypeNameAndDebug(basePointPatchSphericalTensorField, 0);
defineTemplateRunTimeSelectionTable(basePointPatchSphericalTensorField, PointPatch);
defineTemplateRunTimeSelectionTable(basePointPatchSphericalTensorField, patchMapper);
defineTemplateRunTimeSelectionTable(basePointPatchSphericalTensorField, dictionary);

defineNamedTemplateTypeNameAndDebug(basePointPatchSymmTensorField, 0);
defineTemplateRunTimeSelectionTable(basePointPatchSymmTensorField, PointPatch);
defineTemplateRunTimeSelectionTable(basePointPatchSymmTensorField, patchMapper);
defineTemplateRunTimeSelectionTable(basePointPatchSymmTensorField, dictionary);

defineNamedTemplateTypeNameAndDebug(basePointPatchTensorField, 0);
defineTemplateRunTimeSelectionTable(basePointPatchTensorField, PointPatch);
defineTemplateRunTimeSelectionTable(basePointPatchTensorField, patchMapper);
defineTemplateRunTimeSelectionTable(basePointPatchTensorField, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
