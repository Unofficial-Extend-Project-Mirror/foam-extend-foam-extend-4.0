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

#include "polyMesh.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebug(pointScalarField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(pointVectorField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(pointSphericalTensorField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(pointSymmTensorField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(pointSymmTensor4thOrderField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(pointDiagTensorField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(pointTensorField::DimensionedInternalField, 0);

defineTemplateTypeNameAndDebug(pointScalarField, 0);
defineTemplateTypeNameAndDebug(pointVectorField, 0);
defineTemplateTypeNameAndDebug(pointSphericalTensorField, 0);
defineTemplateTypeNameAndDebug(pointSymmTensorField, 0);
defineTemplateTypeNameAndDebug(pointSymmTensor4thOrderField, 0);
defineTemplateTypeNameAndDebug(pointDiagTensorField, 0);
defineTemplateTypeNameAndDebug(pointTensorField, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
