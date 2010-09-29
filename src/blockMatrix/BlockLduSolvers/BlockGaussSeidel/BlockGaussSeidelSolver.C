/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description
    Gauss-Seidel solver for symmetric and asymmetric matrices.  In
    order to improve efficiency, the residual is evaluated after every
    nSweeps sweeps.

\*---------------------------------------------------------------------------*/

#include "BlockGaussSeidelSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from matrix and solver data stream
template<class Type>
Foam::BlockGaussSeidelSolver<Type>::BlockGaussSeidelSolver
(
    const word& fieldName,
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
:
    BlockIterativeSolver<Type>
    (
        fieldName,
        matrix,
        dict
    ),
    gs_(matrix),
    nResSweeps_(readLabel(this->dict().lookup("nResSweeps")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
typename Foam::BlockSolverPerformance<Type>
Foam::BlockGaussSeidelSolver<Type>::solve
(
    Field<Type>& x,
    const Field<Type>& b
)
{
    // Create local references to avoid the spread this-> ugliness
    const BlockLduMatrix<Type>& matrix = this->matrix_;

    // Prepare solver performance
    BlockSolverPerformance<Type> solverPerf
    (
        typeName,
        this->fieldName()
    );

    scalar norm = this->normFactor(x, b);

    Field<Type> wA(x.size());

    // Calculate residual.  Note: sign of residual swapped for efficiency
    matrix.Amul(wA, x);
    wA -= b;

    solverPerf.initialResidual() = gSum(cmptMag(wA))/norm;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // Check convergence, solve if not converged

    if (!solverPerf.checkConvergence(this->tolerance(), this->relTolerance()))
    {
        // Iteration loop

        do
        {
            for (label i = 0; i < nResSweeps_; i++)
            {
                gs_.precondition(x, b);

                solverPerf.nIterations()++;
            }

            // Re-calculate residual.  Note: sign of residual swapped
            // for efficiency
            matrix.Amul(wA, x);
            wA -= b;

            solverPerf.finalResidual() = gSum(cmptMag(wA))/norm;
            solverPerf.nIterations()++;
        } while (!stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
