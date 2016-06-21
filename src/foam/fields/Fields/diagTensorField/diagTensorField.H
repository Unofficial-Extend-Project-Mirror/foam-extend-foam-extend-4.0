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

Typedef
    Foam::diagTensorField

Description
    Specialisation of Field\<T\> for diagTensor.

SourceFiles
    diagTensorField.C

\*---------------------------------------------------------------------------*/

#ifndef diagTensorField_H
#define diagTensorField_H

#include "tensorField.H"
#include "diagTensor.H"

#define TEMPLATE
#include "FieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef Field<diagTensor> diagTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(diagTensor, tensor, diag)
UNARY_FUNCTION(scalar, diagTensor, tr)
UNARY_FUNCTION(sphericalTensor, diagTensor, sph)
UNARY_FUNCTION(scalar, diagTensor, det)
UNARY_FUNCTION(diagTensor, diagTensor, inv)


BINARY_OPERATOR(tensor, diagTensor, tensor, +, add)
BINARY_OPERATOR(tensor, diagTensor, tensor, -, subtract)

BINARY_TYPE_OPERATOR(tensor, diagTensor, tensor, +, add)
BINARY_TYPE_OPERATOR(tensor, diagTensor, tensor, -, subtract)

BINARY_OPERATOR(vector, vector, diagTensor, /, divide)
BINARY_TYPE_OPERATOR(vector, vector, diagTensor, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
