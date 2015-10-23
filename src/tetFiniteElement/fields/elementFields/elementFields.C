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
    Element-based fields for the FEM discretisation.

\*---------------------------------------------------------------------------*/

#include "tetPolyMesh.H"
#include "elementFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebug(elementScalarField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(elementVectorField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(elementSphericalTensorField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(elementSymmTensorField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(elementSymmTensor4thOrderField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(elementDiagTensorField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(elementTensorField::DimensionedInternalField, 0);

defineTemplateTypeNameAndDebug(elementScalarField, 0);
defineTemplateTypeNameAndDebug(elementVectorField, 0);
defineTemplateTypeNameAndDebug(elementSphericalTensorField, 0);
defineTemplateTypeNameAndDebug(elementSymmTensorField, 0);
defineTemplateTypeNameAndDebug(elementSymmTensor4thOrderField, 0);
defineTemplateTypeNameAndDebug(elementDiagTensorField, 0);
defineTemplateTypeNameAndDebug(elementTensorField, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
