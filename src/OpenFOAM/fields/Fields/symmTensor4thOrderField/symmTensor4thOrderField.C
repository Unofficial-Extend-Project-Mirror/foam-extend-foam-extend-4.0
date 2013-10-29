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

#include "symmTensor4thOrderField.H"
#include "transformField.H"
#include "boolList.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

template<>
tmp<Field<symmTensor4thOrder> > transformFieldMask<symmTensor4thOrder>
(
    const symmTensorField& stf
)
{
    tmp<symmTensor4thOrderField> tRes(new symmTensor4thOrderField(stf.size()));
    //symmTensor4thOrderField& res = tRes();
    //TFOR_ALL_F_OP_F(symmTensor4thOrder, res, =, symmTensor, stf)
    return tRes;
}

template<>
tmp<Field<symmTensor4thOrder> > transformFieldMask<symmTensor4thOrder>
(
    const tmp<symmTensorField>& tstf
)
{
    tmp<Field<symmTensor4thOrder> > ret = transformFieldMask<symmTensor4thOrder>(tstf());
    tstf.clear();
    return ret;
}



// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
