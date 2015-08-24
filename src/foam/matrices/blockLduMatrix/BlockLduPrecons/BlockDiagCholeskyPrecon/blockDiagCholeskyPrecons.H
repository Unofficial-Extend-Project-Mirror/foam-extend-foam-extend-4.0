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
    BlockDiagCholeskyPrecon

Description
    Typedefs for Incomplete Cholesky preconditioning with no fill-in,
    using the diagonal-of-diagonal for ILU decomposition.
    Currently broken


Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    blockDiagCholeskyPrecons.C

\*---------------------------------------------------------------------------*/

#ifndef blockDiagCholeskyPrecons_H
#define blockDiagCholeskyPrecons_H

#include "scalarBlockDiagCholeskyPrecon.H"
#include "tensorBlockDiagCholeskyPrecon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef BlockDiagCholeskyPrecon<scalar> blockDiagCholeskyPreconScalar;
typedef BlockDiagCholeskyPrecon<vector> blockDiagCholeskyPreconVector;
typedef BlockDiagCholeskyPrecon<tensor> blockDiagCholeskyPreconTensor;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
