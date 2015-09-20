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
<elementPatchField, elementMesh, tetPolyPatch, DummyMatrix, symmTensor4thOrder>
    baseElementPatchSymmTensor4thOrderField;

typedef PointPatchField
<elementPatchField, elementMesh, tetPolyPatch, DummyMatrix, diagTensor>
    baseElementPatchDiagTensorField;

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

defineNamedTemplateTypeNameAndDebug(baseElementPatchSymmTensor4thOrderField, 0);
defineTemplateRunTimeSelectionTable(baseElementPatchSymmTensor4thOrderField, PointPatch);
defineTemplateRunTimeSelectionTable(baseElementPatchSymmTensor4thOrderField, patchMapper);
defineTemplateRunTimeSelectionTable(baseElementPatchSymmTensor4thOrderField, dictionary);

defineNamedTemplateTypeNameAndDebug(baseElementPatchDiagTensorField, 0);
defineTemplateRunTimeSelectionTable(baseElementPatchDiagTensorField, PointPatch);
defineTemplateRunTimeSelectionTable(baseElementPatchDiagTensorField, patchMapper);
defineTemplateRunTimeSelectionTable(baseElementPatchDiagTensorField, dictionary);

defineNamedTemplateTypeNameAndDebug(baseElementPatchTensorField, 0);
defineTemplateRunTimeSelectionTable(baseElementPatchTensorField, PointPatch);
defineTemplateRunTimeSelectionTable(baseElementPatchTensorField, patchMapper);
defineTemplateRunTimeSelectionTable(baseElementPatchTensorField, dictionary);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
