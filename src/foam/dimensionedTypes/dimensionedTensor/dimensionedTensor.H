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
    Foam::dimensionedTensor

Description
    Dimensioned tensor obtained from generic dimensioned type.

SourceFiles
    dimensionedTensor.C

\*---------------------------------------------------------------------------*/

#ifndef dimensionedTensor_H
#define dimensionedTensor_H

#include "dimensionedVector.H"
#include "dimensionedSymmTensor.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef dimensioned<tensor> dimensionedTensor;


// global functions

dimensionedScalar tr(const dimensionedTensor&);
dimensionedTensor dev(const dimensionedTensor&);
dimensionedTensor dev2(const dimensionedTensor&);
dimensionedScalar det(const dimensionedTensor&);
dimensionedTensor cof(const dimensionedTensor&);
dimensionedTensor inv(const dimensionedTensor&);
dimensionedTensor hinv(const dimensionedTensor&);
dimensionedSymmTensor symm(const dimensionedTensor&);
dimensionedSymmTensor twoSymm(const dimensionedTensor&);
dimensionedTensor skew(const dimensionedTensor&);

dimensionedVector eigenValues(const dimensionedTensor&);
dimensionedTensor eigenVectors(const dimensionedTensor&);

dimensionedVector eigenValues(const dimensionedSymmTensor&);
dimensionedTensor eigenVectors(const dimensionedSymmTensor&);


// global operators

//- Hodge Dual operator (tensor -> vector)
dimensionedVector operator*(const dimensionedTensor&);

//- Hodge Dual operator (vector -> tensor)
dimensionedTensor operator*(const dimensionedVector&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

