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
    Preconditioned Conjugate Gradient solver with run-time selectable
    preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "cgSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cgSolver, 0);

    lduSolver::addsymMatrixConstructorToTable<cgSolver>
        addcgSolverSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
Foam::cgSolver::cgSolver
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& dict
)
:
    lduSolver
    (
        fieldName,
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces,
        dict
    ),
    preconPtr_
    (
        lduPreconditioner::New
        (
            matrix,
            coupleBouCoeffs,
            coupleIntCoeffs,
            interfaces,
            dict
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::lduSolverPerformance Foam::cgSolver::solve
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    // Prepare solver performance
    lduSolverPerformance solverPerf
    (
        preconPtr_->type() + typeName,
        fieldName()
    );

    scalarField wA(x.size());
    scalarField rA(x.size());

    // Calculate initial residual
    matrix_.Amul(wA, x, coupleBouCoeffs_, interfaces_, cmpt);

    // Use rA as scratch space when calculating the normalisation factor
    scalar normFactor = this->normFactor(x, b, wA, rA, cmpt);

    if (lduMatrix::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // Calculate residual
    forAll (rA, i)
    {
        rA[i] = b[i] - wA[i];
    }

    solverPerf.initialResidual() = gSumMag(rA)/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    if (!stop(solverPerf))
    {
        scalar rho = matrix_.great_;
        scalar rhoOld = rho;

        scalar alpha, beta, wApA;

        scalarField pA(x.size(), 0);

        do
        {
            rhoOld = rho;

            // Execute preconditioning
            preconPtr_->precondition(wA, rA, cmpt);

            // Update search directions
            rho = gSumProd(wA, rA);

            beta = rho/rhoOld;

            forAll (pA, i)
            {
                pA[i] = wA[i] + beta*pA[i];
            }


            // Update preconditioned residual
            matrix_.Amul(wA, pA, coupleBouCoeffs_, interfaces_, cmpt);

            wApA = gSumProd(wA, pA);


            // Check for singularity
            if (solverPerf.checkSingularity(mag(wApA)/normFactor))
            {
                break;
            }

            // Update solution and residual
            alpha = rho/wApA;

            forAll (x, i)
            {
                x[i] += alpha*pA[i];
            }

            forAll (rA, i)
            {
                rA[i] -= alpha*wA[i];
            }

            solverPerf.finalResidual() = gSumMag(rA)/normFactor;
            solverPerf.nIterations()++;
        } while (!stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
