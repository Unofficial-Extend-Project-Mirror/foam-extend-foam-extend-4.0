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
    areaFields

Description

SourceFiles
    areaFields.C


\*---------------------------------------------------------------------------*/

#ifndef areaFieldsFwd_H
#define areaFieldsFwd_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class areaMesh;

template<class Type>
class faPatchField;

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField;

typedef GeometricField<scalar, faPatchField, areaMesh> areaScalarField;
typedef GeometricField<vector, faPatchField, areaMesh> areaVectorField;
typedef GeometricField<sphericalTensor, faPatchField, areaMesh>
    areaSphericalTensorField;
typedef GeometricField<symmTensor, faPatchField, areaMesh> areaSymmTensorField;
typedef GeometricField<tensor, faPatchField, areaMesh> areaTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
