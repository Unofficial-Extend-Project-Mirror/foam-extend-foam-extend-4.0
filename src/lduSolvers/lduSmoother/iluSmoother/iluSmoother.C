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
    iluSmoother

Description
    Symmetric Gauss-Seidel smoother

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "iluSmoother.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(iluSmoother, 0);

    lduSmoother::addsymMatrixConstructorToTable<iluSmoother>
        addiluSmootherSymMatrixConstructorToTable_;

    lduSmoother::addasymMatrixConstructorToTable<iluSmoother>
        addiluSmootherAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::iluSmoother::iluSmoother
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    lduSmoother
    (
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces
    ),
    precon_
    (
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces
    ),
    xCorr_(matrix.lduAddr().size()),
    residual_(matrix.lduAddr().size())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::iluSmoother::smooth
(
    scalarField& x,
    const scalarField& b,
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
            coupleBouCoeffs_,
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
