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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Algebraic Multigrid solver with run-time selection of policy and cycle

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "amgSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(amgSolver, 0);

    lduSolver::addsymMatrixConstructorToTable<amgSolver>
        addamgSolverSymMatrixConstructorToTable_;

    lduSolver::addasymMatrixConstructorToTable<amgSolver>
        addamgSolverAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
Foam::amgSolver::amgSolver
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
    amg_
    (
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces,
        dict
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::lduSolverPerformance Foam::amgSolver::solve
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    // Prepare solver performance
    lduSolverPerformance solverPerf(typeName, fieldName());

    scalar norm = this->normFactor(x, b, cmpt);

    // Calculate initial residual
    solverPerf.initialResidual() = gSumMag(amg_.residual(x, b, cmpt))/norm;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    if (!stop(solverPerf))
    {
        do
        {
            amg_.cycle(x, b, cmpt);

            // Re-calculate residual
            solverPerf.finalResidual() =
                gSumMag(amg_.residual(x, b, cmpt))/norm;

            solverPerf.nIterations()++;
        } while (!stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
