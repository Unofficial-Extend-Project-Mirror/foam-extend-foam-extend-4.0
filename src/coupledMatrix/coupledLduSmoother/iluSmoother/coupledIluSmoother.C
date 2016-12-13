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
    Symmetric Gauss-Seidel smoother for coupled lduMatrices.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "coupledIluSmoother.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledIluSmoother, 0);

    addToRunTimeSelectionTable
    (
        coupledLduSmoother,
        coupledIluSmoother,
        word
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::coupledIluSmoother::coupledIluSmoother
(
    const coupledLduMatrix& matrix,
    const PtrList<FieldField<Field, scalar> >& bouCoeffs,
    const PtrList<FieldField<Field, scalar> >& intCoeffs,
    const lduInterfaceFieldPtrsListList& interfaces
)
:
    coupledLduSmoother
    (
        matrix,
        bouCoeffs,
        intCoeffs,
        interfaces
    ),
    precon_
    (
        matrix,
        bouCoeffs,
        intCoeffs,
        interfaces
    ),
    xCorr_(matrix.size()),
    residual_(matrix.size())
{
    forAll (xCorr_, rowI)
    {
        xCorr_.set(rowI, new scalarField(matrix[rowI].lduAddr().size(), 0));
        residual_.set(rowI, new scalarField(matrix[rowI].lduAddr().size(), 0));
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::coupledIluSmoother::smooth
(
    FieldField<Field, scalar>& x,
    const FieldField<Field, scalar>& b,
    const direction cmpt,
    const label nSweeps
) const
{
    for (label sweep = 0; sweep < nSweeps; sweep++)
    {
        // Calculate residual
        matrix_.Amul
        (
            residual_,
            x,
            bouCoeffs_,
            interfaces_,
            cmpt
        );

        // residual = b - Ax
        forAll (b, i)
        {
            residual_[i] = b[i] - residual_[i];
        }

        precon_.precondition(xCorr_, residual_, cmpt);

        // Add correction to x
        x += xCorr_;
    }
}


// ************************************************************************* //
