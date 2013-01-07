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

Description
    Preconditioned Bi-Conjugate Gradient stabilised solver with run-time
    selectable preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "bicgStabSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bicgStabSolver, 0);

    //HJ, temporary
    lduSolver::addsymMatrixConstructorToTable<bicgStabSolver>
        addbicgStabSolverSymMatrixConstructorToTable_;

    lduSolver::addasymMatrixConstructorToTable<bicgStabSolver>
        addbicgStabSolverAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
Foam::bicgStabSolver::bicgStabSolver
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

Foam::lduSolverPerformance Foam::bicgStabSolver::solve
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    // Prepare solver performance
    lduSolverPerformance solverPerf(typeName, fieldName());

    scalarField p(x.size());
    scalarField r(x.size());

    // Calculate initial residual
    matrix_.Amul(p, x, coupleBouCoeffs_, interfaces_, cmpt);

    scalar normFactor = this->normFactor(x, b, p, r, cmpt);

    if (lduMatrix::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // Calculate residual
    forAll (r, i)
    {
        r[i] = b[i] - p[i];
    }

    solverPerf.initialResidual() = gSumMag(r)/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    if (!stop(solverPerf))
    {
        scalar rho = matrix_.great_;
        scalar rhoOld = rho;

        scalar alpha = 0;
        scalar omega = matrix_.great_;
        scalar beta;

        p = 0;
        scalarField ph(x.size(), 0);
        scalarField v(x.size(), 0);
        scalarField s(x.size(), 0);
        scalarField sh(x.size(), 0);
        scalarField t(x.size(), 0);

        // Calculate transpose residual
        scalarField rw(r);

        do
        {
            rhoOld = rho;

            // Update search directions
            rho = gSumProd(rw, r);

            beta = rho/rhoOld*(alpha/omega);

            // Restart if breakdown occurs
            if (rho == 0)
            {
                rw = r;
                rho = gSumProd(rw, r);

                alpha = 0;
                omega = 0;
                beta = 0;
            }

            forAll (p, i)
            {
                p[i] = r[i] + beta*p[i] - beta*omega*v[i];
            }

            // Execute preconditioning
            preconPtr_->precondition(ph, p, cmpt);
            matrix_.Amul(v, ph, coupleBouCoeffs_, interfaces_, cmpt);
            alpha = rho/gSumProd(rw, v);

            forAll (s, i)
            {
                s[i] = r[i] - alpha*v[i];
            }

            // Execute preconditioning
            // Bug fix, Alexander Monakov, 11/Jul/2012
            preconPtr_->precondition(sh, s, cmpt);
            matrix_.Amul(t, sh, coupleBouCoeffs_, interfaces_, cmpt);
            omega = gSumProd(t, s)/gSumProd(t, t);

            // Update solution and residual
            forAll (x, i)
            {
                x[i] = x[i] + alpha*ph[i] + omega*sh[i];
            }

            forAll (r, i)
            {
                r[i] = s[i] - omega*t[i];
            }

            solverPerf.finalResidual() = gSumMag(r)/normFactor;
            solverPerf.nIterations()++;
        } while (!stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
