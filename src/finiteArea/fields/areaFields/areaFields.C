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
#include "areaFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebug(areaScalarField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(areaVectorField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(areaSphericalTensorField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(areaSymmTensorField::DimensionedInternalField, 0);
defineTemplateTypeNameAndDebug(areaTensorField::DimensionedInternalField, 0);

defineTemplateTypeNameAndDebug(areaScalarField, 0);
defineTemplateTypeNameAndDebug(areaVectorField, 0);
defineTemplateTypeNameAndDebug(areaSphericalTensorField, 0);
defineTemplateTypeNameAndDebug(areaSymmTensorField, 0);
defineTemplateTypeNameAndDebug(areaTensorField, 0);

template<>
tmp<GeometricField<scalar, faPatchField, areaMesh> >
GeometricField<scalar, faPatchField, areaMesh>::component
(
    const direction
) const
{
    return *this;
}

template<>
void GeometricField<scalar, faPatchField, areaMesh>::replace
(
    const direction,
    const GeometricField<scalar, faPatchField, areaMesh>& gsf
)
{
    *this == gsf;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

