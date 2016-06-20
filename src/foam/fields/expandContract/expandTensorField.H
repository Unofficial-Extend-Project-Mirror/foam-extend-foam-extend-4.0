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
    Global functions for expansion and contraction of tensor field
    to diagonal type

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    expandTensorField.C

\*---------------------------------------------------------------------------*/

#ifndef expandTensorField_H
#define expandTensorField_H

#include "tensorField.H"

#define TEMPLATE
#include "FieldFunctionsM.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void contractScalar(Field<scalar>& res, const UList<tensor>& f);
void contractLinear(Field<vector>& res, const UList<tensor>& f);

void expandScalar(Field<vector>& res, const UList<scalar>& f);

void expandScalar(Field<tensor>& res, const UList<scalar>& f);
void expandLinear(Field<tensor>& res, const UList<vector>& f);

void sumToDiag(Field<vector>& res, const UList<tensor>& f);
void sumMagToDiag(Field<vector>& res, const UList<tensor>& f);

} // End namespace Foam

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
