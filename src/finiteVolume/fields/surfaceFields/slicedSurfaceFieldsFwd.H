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

InClass
    Foam::slicedSurfaceFields

\*---------------------------------------------------------------------------*/

#ifndef slicedSurfaceFieldsFwd_H
#define slicedSurfaceFieldsFwd_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class surfaceMesh;

template<class Type>
class fvsPatchField;

template<class Type>
class slicedFvsPatchField;

template
<
    class Type,
    template<class> class PatchField,
    template<class> class SlicedPatchField,
    class GeoMesh
>
class SlicedGeometricField;

typedef
SlicedGeometricField<scalar, fvsPatchField, slicedFvsPatchField, surfaceMesh>
    slicedSurfaceScalarField;

typedef
SlicedGeometricField<vector, fvsPatchField, slicedFvsPatchField, surfaceMesh>
    slicedSurfaceVectorField;

typedef
SlicedGeometricField<sphericalTensor, fvsPatchField, slicedFvsPatchField, surfaceMesh>
    slicedSurfaceSphericalTensorField;

typedef
SlicedGeometricField<symmTensor, fvsPatchField, slicedFvsPatchField, surfaceMesh>
    slicedSurfaceSymmTensorField;

typedef
SlicedGeometricField<symmTensor4thOrder, fvsPatchField, slicedFvsPatchField, surfaceMesh>
    slicedSurfaceSymmTensor4thOrderField;

typedef
SlicedGeometricField<diagTensor, fvsPatchField, slicedFvsPatchField, surfaceMesh>
    slicedSurfaceDiagTensorField;

typedef
SlicedGeometricField<tensor, fvsPatchField, slicedFvsPatchField, surfaceMesh>
    slicedSurfaceTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
