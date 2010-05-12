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
defineTemplateTypeNameAndDebug(tetPointTensorField::DimensionedInternalField, 0);

defineTemplateTypeNameAndDebug(tetPointScalarField, 0);
defineTemplateTypeNameAndDebug(tetPointVectorField, 0);
defineTemplateTypeNameAndDebug(tetPointSphericalTensorField, 0);
defineTemplateTypeNameAndDebug(tetPointSymmTensorField, 0);
defineTemplateTypeNameAndDebug(tetPointTensorField, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
