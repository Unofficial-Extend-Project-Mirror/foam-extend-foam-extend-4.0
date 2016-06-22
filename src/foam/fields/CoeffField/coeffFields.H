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

Class
    coeffFields

Description
    Typedefs for coefficient field

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef coeffFields_H
#define coeffFields_H

#include "scalarCoeffField.H"
#include "sphericalTensorCoeffField.H"
#include "diagTensorCoeffField.H"
#include "symmTensorCoeffField.H"
#include "tensorCoeffField.H"
#include "symmTensor4thOrderCoeffField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef CoeffField<scalar> scalarCoeffField;
typedef CoeffField<vector> vectorCoeffField;
typedef CoeffField<sphericalTensor> sphericalTensorCoeffField;
typedef CoeffField<symmTensor> diagTensorCoeffField;
typedef CoeffField<symmTensor> symmTensorCoeffField;
typedef CoeffField<tensor> tensorCoeffField;
typedef CoeffField<symmTensor4thOrder> symmTensor4thOrderCoeffField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
