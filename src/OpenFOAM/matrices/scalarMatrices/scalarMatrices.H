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

Class
    scalarMatrices

Description
    Scalar matrices

SourceFiles
    scalarMatrices.C
    scalarMatricesTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef scalarMatrices_H
#define scalarMatrices_H

#include "RectangularMatrix.H"
#include "SquareMatrix.H"
#include "DiagonalMatrix.H"

#include "scalarSquareMatrix.H"
#include "scalarField.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef RectangularMatrix<scalar> scalarRectangularMatrix;
typedef DiagonalMatrix<scalar> scalarDiagonalMatrix;

void multiply
(
    scalarRectangularMatrix& answer,         // value changed in return
    const scalarRectangularMatrix& A,
    const scalarRectangularMatrix& B
);

void multiply
(
    scalarRectangularMatrix& answer,         // value changed in return
    const scalarRectangularMatrix& A,
    const scalarRectangularMatrix& B,
    const scalarRectangularMatrix& C
);

void multiply
(
    scalarRectangularMatrix& answer,         // value changed in return
    const scalarRectangularMatrix& A,
    const scalarDiagonalMatrix& B,
    const scalarRectangularMatrix& C
);

//- Return the inverse of matrix A using SVD
scalarRectangularMatrix SVDinv
(
    const scalarRectangularMatrix& A,
    scalar minCondition = 0
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
