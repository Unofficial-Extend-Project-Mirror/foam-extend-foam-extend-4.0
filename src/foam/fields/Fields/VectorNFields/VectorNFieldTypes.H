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

Description
    Include the header files for all the primitive types that Fields are
    instantiated for.

\*---------------------------------------------------------------------------*/

#ifndef VectorNFieldTypes_H
#define VectorNFieldTypes_H

#include "vector2.H"
#include "vector3.H"
#include "vector4.H"
#include "vector6.H"
#include "vector8.H"

#include "tensor2.H"
#include "tensor3.H"
#include "tensor4.H"
#include "tensor6.H"
#include "tensor8.H"

#include "diagTensor2.H"
#include "diagTensor3.H"
#include "diagTensor4.H"
#include "diagTensor6.H"
#include "diagTensor8.H"

#define forAllVectorNTypes(m, args...)              \
    m(vector2, Vector2, args)                       \
    m(vector3, Vector3, args)                       \
    m(vector4, Vector4, args)                       \
    m(vector6, Vector6, args)                       \
    m(vector8, Vector8, args)

#define forAllTensorNTypes(m, args...)              \
    m(tensor2, Tensor2, args)                       \
    m(tensor3, Tensor3, args)                       \
    m(tensor4, Tensor4, args)                       \
    m(tensor6, Tensor6, args)                       \
    m(tensor8, Tensor8, args)

#define forAllDiagTensorNTypes(m, args...)          \
    m(diagTensor2, DiagTensor2, args)               \
    m(diagTensor3, DiagTensor3, args)               \
    m(diagTensor4, DiagTensor4, args)               \
    m(diagTensor6, DiagTensor6, args)               \
    m(diagTensor8, DiagTensor8, args)

#define forAllSphericalTensorNTypes(m, args...)     \
    m(sphericalTensor2, SphericalTensor2, args)     \
    m(sphericalTensor3, SphericalTensor3, args)     \
    m(sphericalTensor4, SphericalTensor4, args)     \
    m(sphericalTensor6, SphericalTensor6, args)     \
    m(sphericalTensor8, SphericalTensor8, args)

#define forAllVectorTensorNTypes(m, args...)                            \
    m(tensor2, diagTensor2, sphericalTensor2, vector2, scalar, args)    \
    m(tensor3, diagTensor3, sphericalTensor3, vector3, scalar, args)    \
    m(tensor4, diagTensor4, sphericalTensor4, vector4, scalar, args)    \
    m(tensor6, diagTensor6, sphericalTensor6, vector6, scalar, args)    \
    m(tensor8, diagTensor8, sphericalTensor8, vector8, scalar, args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
