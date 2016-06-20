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
<pointPatchField, pointMesh, pointPatch, DummyMatrix, symmTensor4thOrder>
    basePointPatchSymmTensor4thOrderField;

typedef PointPatchField
<pointPatchField, pointMesh, pointPatch, DummyMatrix, diagTensor>
    basePointPatchDiagTensorField;

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

defineNamedTemplateTypeNameAndDebug(basePointPatchSymmTensor4thOrderField, 0);
defineTemplateRunTimeSelectionTable(basePointPatchSymmTensor4thOrderField, PointPatch);
defineTemplateRunTimeSelectionTable(basePointPatchSymmTensor4thOrderField, patchMapper);
defineTemplateRunTimeSelectionTable(basePointPatchSymmTensor4thOrderField, dictionary);

defineNamedTemplateTypeNameAndDebug(basePointPatchDiagTensorField, 0);
defineTemplateRunTimeSelectionTable(basePointPatchDiagTensorField, PointPatch);
defineTemplateRunTimeSelectionTable(basePointPatchDiagTensorField, patchMapper);
defineTemplateRunTimeSelectionTable(basePointPatchDiagTensorField, dictionary);

defineNamedTemplateTypeNameAndDebug(basePointPatchTensorField, 0);
defineTemplateRunTimeSelectionTable(basePointPatchTensorField, PointPatch);
defineTemplateRunTimeSelectionTable(basePointPatchTensorField, patchMapper);
defineTemplateRunTimeSelectionTable(basePointPatchTensorField, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
