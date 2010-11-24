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
    Specialisation of Field<T> for given FadOne<nVars>.

\*---------------------------------------------------------------------------*/

#include "fadOneFields.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

BINARY_TYPE_OPERATOR(fadScalar, fadScalar, fadScalar, +, add)
BINARY_TYPE_OPERATOR(fadScalar, fadScalar, fadScalar, -, subtract)

BINARY_OPERATOR(fadScalar, fadScalar, fadScalar, *, multiply)
BINARY_OPERATOR(fadScalar, fadScalar, fadScalar, /, divide)

BINARY_TYPE_OPERATOR_SF(fadScalar, fadScalar, fadScalar, /, divide)

BINARY_FUNCTION(fadScalar, fadScalar, fadScalar, pow)
BINARY_TYPE_FUNCTION(fadScalar, fadScalar, fadScalar, pow)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(fadScalar, fadScalar, pos)
UNARY_FUNCTION(fadScalar, fadScalar, neg)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<fadScalarField>
fadScalarField::component(const direction) const
{
    return *this;
}

template<>
void component
(
    fadScalarField& lf,
    const UList<fadScalar>& f,
    const direction
)
{
    lf = f;
}

template<>
void fadScalarField::replace
(
    const direction,
    const UList<fadScalar>& lf
)
{
    *this = lf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
