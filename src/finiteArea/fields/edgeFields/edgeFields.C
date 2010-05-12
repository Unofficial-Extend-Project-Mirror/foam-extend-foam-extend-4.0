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

#include "faMesh.H"
#include "edgeFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebug(edgeScalarField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(edgeVectorField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(edgeSphericalTensorField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(edgeSymmTensorField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(edgeTensorField::DimensionedInternalField, 0);

defineTemplateTypeNameAndDebug(edgeScalarField, 0);
defineTemplateTypeNameAndDebug(edgeVectorField, 0);
defineTemplateTypeNameAndDebug(edgeSphericalTensorField, 0);
defineTemplateTypeNameAndDebug(edgeSymmTensorField, 0);
defineTemplateTypeNameAndDebug(edgeTensorField, 0);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

