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
    labelField

Description
    Specialisation of Field<T> for given FadOne<nVars>.

SourceFiles
    fadOneFields.C

\*---------------------------------------------------------------------------*/

#ifndef fadOneFields_H
#define fadOneFields_H

#include "FadOneField.H"

#define TEMPLATE
#include "FieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef FadOne<2> fadScalar;
typedef Field<fadScalar> fadScalarField;

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
fadScalarField::component
(
    const direction
) const;


template<>
void component
(
    fadScalarField& lf,
    const UList<fadScalar>& f,
    const direction
);

template<>
void fadScalarField::replace
(
    const direction,
    const UList<fadScalar>& lf
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
