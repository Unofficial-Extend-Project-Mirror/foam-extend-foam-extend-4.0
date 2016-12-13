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
    Foam::slicedVolFields

\*---------------------------------------------------------------------------*/

#ifndef slicedVolFieldsFwd_H
#define slicedVolFieldsFwd_H

#include "fvPatchField.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class volMesh;

template<class Type>
class fvPatchField;

template<class Type>
class slicedFvPatchField;

template
<
    class Type,
    template<class> class PatchField,
    template<class> class SlicedPatchField,
    class GeoMesh
>
class SlicedGeometricField;

typedef
SlicedGeometricField<scalar, fvPatchField, slicedFvPatchField, volMesh>
    slicedVolScalarField;

typedef
SlicedGeometricField<vector, fvPatchField, slicedFvPatchField, volMesh>
    slicedVolVectorField;

typedef
SlicedGeometricField<sphericalTensor, fvPatchField, slicedFvPatchField, volMesh>
    slicedVolSphericalTensorField;

typedef
SlicedGeometricField<symmTensor, fvPatchField, slicedFvPatchField, volMesh>
    slicedVolSymmTensorField;

typedef
SlicedGeometricField<symmTensor4thOrder, fvPatchField, slicedFvPatchField, volMesh>
    slicedVolSymmTensor4thOrderField;

typedef
SlicedGeometricField<diagTensor, fvPatchField, slicedFvPatchField, volMesh>
    slicedVolDiagTensorField;

typedef
SlicedGeometricField<tensor, fvPatchField, slicedFvPatchField, volMesh>
    slicedVolTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
