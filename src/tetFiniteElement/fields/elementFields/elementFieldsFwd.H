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

Class
    elementFields

Description
    Element-based fields for the FEM discretisation.

SourceFiles
    elementFields.C


\*---------------------------------------------------------------------------*/

#ifndef elementFieldsFwd_H
#define elementFieldsFwd_H

#include "fieldTypes.H"
#include "elementMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
class elementPatchField;

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField;

typedef
GeometricField<scalar, elementPatchField, elementMesh> elementScalarField;

typedef
GeometricField<vector, elementPatchField, elementMesh> elementVectorField;

typedef
GeometricField<sphericalTensor, elementPatchField, elementMesh>
    elementSphericalTensorField;

typedef
GeometricField<symmTensor, elementPatchField, elementMesh>
    elementSymmTensorField;

typedef
GeometricField<tensor, elementPatchField, elementMesh> elementTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
