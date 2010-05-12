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
    edgeFields

Description

SourceFiles
    edgeFields.C


\*---------------------------------------------------------------------------*/

#ifndef edgeFieldsFwd_H
#define edgeFieldsFwd_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class edgeMesh;

template<class Type>
class faePatchField;

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField;

typedef GeometricField<scalar, faePatchField, edgeMesh> edgeScalarField;
typedef GeometricField<vector, faePatchField, edgeMesh> edgeVectorField;
typedef GeometricField<sphericalTensor, faePatchField, edgeMesh>
    edgeSphericalTensorField;
typedef GeometricField<symmTensor, faePatchField, edgeMesh> edgeSymmTensorField;
typedef GeometricField<tensor, faePatchField, edgeMesh> edgeTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
