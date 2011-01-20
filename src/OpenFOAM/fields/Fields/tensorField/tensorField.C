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

\*---------------------------------------------------------------------------*/

#include "tensorField.H"
#include "transformField.H"
#include "boolList.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, tensor, tr)
UNARY_FUNCTION(sphericalTensor, tensor, sph)
UNARY_FUNCTION(symmTensor, tensor, symm)
UNARY_FUNCTION(symmTensor, tensor, twoSymm)
UNARY_FUNCTION(tensor, tensor, skew)
UNARY_FUNCTION(tensor, tensor, dev)
UNARY_FUNCTION(tensor, tensor, dev2)
UNARY_FUNCTION(scalar, tensor, det)
UNARY_FUNCTION(tensor, tensor, cof)

UNARY_FUNCTION(tensor, tensor, hinv)


// This is a nasty hack for flat geometries.  For serious SVD, please use hinv
// (Singular value decomposition and Householder transfromations)
// HJ, 24/Oct/2009
void inv(Field<tensor>& tf, const UList<tensor>& tf1)
{
    if (tf.empty())
    {
        return;
    }

    scalar scale = magSqr(tf1[0]);

    // Fixed terrible hack.  HJ, 20/Jan/2011
    boolList removeCmpts(3);
    removeCmpts[0] = magSqr(tf1[0].xx())/scale < SMALL;
    removeCmpts[1] = magSqr(tf1[0].yy())/scale < SMALL;
    removeCmpts[2] = magSqr(tf1[0].zz())/scale < SMALL;


    if (removeCmpts[0] || removeCmpts[1] || removeCmpts[2])
    {
        tensorField tf1Plus(tf1);

        if (removeCmpts[0])
        {
            tf1Plus += tensor(1,0,0,0,0,0,0,0,0);
        }

        if (removeCmpts[1])
        {
            tf1Plus += tensor(0,0,0,0,1,0,0,0,0);
        }

        if (removeCmpts[2])
        {
            tf1Plus += tensor(0,0,0,0,0,0,0,0,1);
        }

        TFOR_ALL_F_OP_FUNC_F(tensor, tf, =, inv, tensor, tf1Plus)

        if (removeCmpts[0])
        {
            tf -= tensor(1,0,0,0,0,0,0,0,0);
        }

        if (removeCmpts[1])
        {
            tf -= tensor(0,0,0,0,1,0,0,0,0);
        }

        if (removeCmpts[2])
        {
            tf -= tensor(0,0,0,0,0,0,0,0,1);
        }
    }
    else
    {
        TFOR_ALL_F_OP_FUNC_F(tensor, tf, =, inv, tensor, tf1)
    }
}


tmp<tensorField> inv(const UList<tensor>& tf)
{
    tmp<tensorField> result(new tensorField(tf.size()));
    inv(result(), tf);
    return result;
}


tmp<tensorField> inv(const tmp<tensorField>& tf)
{
    tmp<tensorField> tRes = reuseTmp<tensor, tensor>::New(tf);
    inv(tRes(), tf());
    reuseTmp<tensor, tensor>::clear(tf);
    return tRes;
}


UNARY_FUNCTION(vector, tensor, eigenValues)
UNARY_FUNCTION(tensor, tensor, eigenVectors)

UNARY_FUNCTION(vector, symmTensor, eigenValues)
UNARY_FUNCTION(tensor, symmTensor, eigenVectors)


template<>
tmp<Field<tensor> > transformFieldMask<tensor>
(
    const symmTensorField& stf
)
{
    tmp<tensorField> tRes(new tensorField(stf.size()));
    tensorField& res = tRes();
    TFOR_ALL_F_OP_F(tensor, res, =, symmTensor, stf)
    return tRes;
}

template<>
tmp<Field<tensor> > transformFieldMask<tensor>
(
    const tmp<symmTensorField>& tstf
)
{
    tmp<Field<tensor> > ret = transformFieldMask<tensor>(tstf());
    tstf.clear();
    return ret;
}


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, tensor, *, hdual)
UNARY_OPERATOR(tensor, vector, *, hdual)

BINARY_OPERATOR(vector, vector, tensor, /, divide)
BINARY_TYPE_OPERATOR(vector, vector, tensor, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
