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
    dimensionedVectorTensorN

Description
    Dimensioned VectorN and TensorN obtained from generic dimensioned type.

Author
    Ivor Clifford

SourceFiles
    dimensionedVectorTensorN.C

\*---------------------------------------------------------------------------*/

#ifndef dimensionedVectorTensorN_H
#define dimensionedVectorTensorN_H

#include "dimensionedType.H"
#include "VectorNFieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

#define doMakeTypedef(type, Type, args...)      \
    typedef dimensioned<type> dimensioned##Type;

forAllVectorNTypes(doMakeTypedef)

forAllTensorNTypes(doMakeTypedef)

forAllDiagTensorNTypes(doMakeTypedef)

forAllSphericalTensorNTypes(doMakeTypedef)

#undef doMakeTypedef


// friend functions

dimensionSet diag(const dimensionSet& ds);
dimensionSet cmptSum(const dimensionSet& ds);
dimensionSet contractLinear(const dimensionSet& ds);
dimensionSet contractScalar(const dimensionSet& ds);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define dimensionedType_Funs(tensorType, diagTensorType,        \
    sphericalTensorType, vectorType, cmptType, args...)         \
                                                                \
tensorType inv(const tensorType&);                              \
diagTensorType inv(const diagTensorType&);                      \
sphericalTensorType inv(const sphericalTensorType&);            \
                                                                \
tensorType inv(const tensorType&);                              \
diagTensorType inv(const diagTensorType&);                      \
sphericalTensorType inv(const sphericalTensorType&);            \
                                                                \
cmptType cmptSum(const vectorType&);

#define dimensionedVectorN_Funs(length)                         \
dimensionedType_Funs                                            \
(                                                               \
    dimensionedScalar,                                          \
    dimensionedVector##length,                                  \
    dimensionedTensor##length,                                  \
    dimensionedDiagTensor##length,                              \
    dimensionedSphericalTensor##length                          \
)

forAllVectorTensorNTypes(dimensionedType_Funs)

#undef dimensionedVectorN_Funs
#undef dimensionedType_Funs

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dimensionedVectorTensorNI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

