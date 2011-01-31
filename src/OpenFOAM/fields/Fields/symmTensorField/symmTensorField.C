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

#include "symmTensorField.H"
#include "transformField.H"
#include "boolList.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

UNARY_FUNCTION(symmTensor, vector, sqr)

UNARY_FUNCTION(scalar, symmTensor, tr)
UNARY_FUNCTION(sphericalTensor, symmTensor, sph)
UNARY_FUNCTION(symmTensor, symmTensor, symm)
UNARY_FUNCTION(symmTensor, symmTensor, twoSymm)
UNARY_FUNCTION(symmTensor, symmTensor, dev)
UNARY_FUNCTION(symmTensor, symmTensor, dev2)
UNARY_FUNCTION(scalar, symmTensor, det)
UNARY_FUNCTION(symmTensor, symmTensor, cof)

UNARY_FUNCTION(symmTensor, symmTensor, hinv)

// This is a nasty hack for flat geometries.  For serious SVD, please use hinv
// (Singular value decomposition and Householder transfromations)
// HJ, 24/Oct/2009
void inv(Field<symmTensor>& tf, const UList<symmTensor>& tf1)
{
    if (tf.empty())
    {
        return;
    }

    scalar scale = magSqr(tf1[0]);

    // Fixed terrible hack.  HJ, 20/Jan/2011
    boolList removeCmpts(3);
    removeCmpts[0] =  magSqr(tf1[0].xx())/scale < SMALL;
    removeCmpts[1] =  magSqr(tf1[0].yy())/scale < SMALL;
    removeCmpts[2] =  magSqr(tf1[0].zz())/scale < SMALL;

    if (removeCmpts[0] || removeCmpts[1] || removeCmpts[2])
    {
        symmTensorField tf1Plus(tf1);

        if (removeCmpts[0])
        {
            tf1Plus += symmTensor(1,0,0,0,0,0);
        }

        if (removeCmpts[1])
        {
            tf1Plus += symmTensor(0,0,0,1,0,0);
        }

        if (removeCmpts[2])
        {
            tf1Plus += symmTensor(0,0,0,0,0,1);
        }

        TFOR_ALL_F_OP_FUNC_F(symmTensor, tf, =, inv, symmTensor, tf1Plus)

        if (removeCmpts[0])
        {
            tf -= symmTensor(1,0,0,0,0,0);
        }

        if (removeCmpts[1])
        {
            tf -= symmTensor(0,0,0,1,0,0);
        }

        if (removeCmpts[2])
        {
            tf -= symmTensor(0,0,0,0,0,1);
        }
    }
    else
    {
        TFOR_ALL_F_OP_FUNC_F(symmTensor, tf, =, inv, symmTensor, tf1)
    }
}

tmp<symmTensorField> inv(const UList<symmTensor>& tf)
{
    tmp<symmTensorField> result(new symmTensorField(tf.size()));
    inv(result(), tf);
    return result;
}


tmp<symmTensorField> inv(const tmp<symmTensorField>& tf)
{
    tmp<symmTensorField> tRes = reuseTmp<symmTensor, symmTensor>::New(tf);
    inv(tRes(), tf());
    reuseTmp<symmTensor, symmTensor>::clear(tf);
    return tRes;
}


template<>
tmp<Field<symmTensor> > transformFieldMask<symmTensor>
(
    const tensorField& tf
)
{
    return symm(tf);
}

template<>
tmp<Field<symmTensor> > transformFieldMask<symmTensor>
(
    const tmp<tensorField>& ttf
)
{
    tmp<Field<symmTensor> > ret = transformFieldMask<symmTensor>(ttf());
    ttf.clear();
    return ret;
}


template<>
tmp<Field<symmTensor> > transformFieldMask<symmTensor>
(
    const symmTensorField& stf
)
{
    return stf;
}

template<>
tmp<Field<symmTensor> > transformFieldMask<symmTensor>
(
    const tmp<symmTensorField>& tstf
)
{
    return tstf;
}


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, symmTensor, *, hdual)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
