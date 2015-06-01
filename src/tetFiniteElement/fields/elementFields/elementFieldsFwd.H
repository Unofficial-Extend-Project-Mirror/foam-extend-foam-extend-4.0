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
GeometricField<symmTensor4thOrder, elementPatchField, elementMesh>
    elementSymmTensor4thOrderField;

typedef
GeometricField<diagTensor, elementPatchField, elementMesh>
    elementDiagTensorField;

typedef
GeometricField<tensor, elementPatchField, elementMesh> elementTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
