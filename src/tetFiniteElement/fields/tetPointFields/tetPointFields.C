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
    Point-based fields for the FEM discretisation.

\*---------------------------------------------------------------------------*/

#include "tetPolyMesh.H"
#include "tetPointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebug(tetPointScalarField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(tetPointVectorField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(tetPointSphericalTensorField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(tetPointSymmTensorField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(tetPointSymmTensor4thOrderField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(tetPointDiagTensorField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(tetPointTensorField::DimensionedInternalField, 0);

defineTemplateTypeNameAndDebug(tetPointScalarField, 0);
defineTemplateTypeNameAndDebug(tetPointVectorField, 0);
defineTemplateTypeNameAndDebug(tetPointSphericalTensorField, 0);
defineTemplateTypeNameAndDebug(tetPointSymmTensorField, 0);
defineTemplateTypeNameAndDebug(tetPointSymmTensor4thOrderField, 0);
defineTemplateTypeNameAndDebug(tetPointDiagTensorField, 0);
defineTemplateTypeNameAndDebug(tetPointTensorField, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
