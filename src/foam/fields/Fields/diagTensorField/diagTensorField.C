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
    Specialisation of Field\<T\> for diagTensor.

\*---------------------------------------------------------------------------*/

#include "diagTensorField.H"
#include "transformField.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

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


template<>
tmp<Field<diagTensor> > transformFieldMask<diagTensor>
(
    const tensorField& tf
)
{
    tmp<Field<diagTensor> > ret(new Field<diagTensor>(tf.size()));

    ret().component(diagTensor::XX) = tf.component(tensor::XX);
    ret().component(diagTensor::YY) = tf.component(tensor::YY);
    ret().component(diagTensor::ZZ) = tf.component(tensor::ZZ);

    return ret;
}

template<>
tmp<Field<diagTensor> > transformFieldMask<diagTensor>
(
    const tmp<tensorField>& ttf
)
{
    tmp<Field<diagTensor> > ret =
        transformFieldMask<diagTensor>(ttf());
    ttf.clear();
    return ret;
}


template<>
tmp<Field<diagTensor> > transformFieldMask<diagTensor>
(
    const symmTensorField& stf
)
{
    tmp<Field<diagTensor> > ret(new Field<diagTensor>(stf.size()));

    ret().component(diagTensor::XX) = stf.component(symmTensor::XX);
    ret().component(diagTensor::YY) = stf.component(symmTensor::YY);
    ret().component(diagTensor::ZZ) = stf.component(symmTensor::ZZ);

    return ret;
}

template<>
tmp<Field<diagTensor> > transformFieldMask<diagTensor>
(
    const tmp<symmTensorField>& tstf
)
{
    tmp<Field<diagTensor> > ret =
        transformFieldMask<diagTensor>(tstf());
    tstf.clear();
    return ret;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
