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
    Forward Projection Euler Algebraic Multigrid solver with run-time selection
    of policy and cycle

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "fpeAmgSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(fpeAmgSolver, 0);

    lduSolver::addsymMatrixConstructorToTable<fpeAmgSolver>
        addfpeAmgSolverSymMatrixConstructorToTable_;

    lduSolver::addasymMatrixConstructorToTable<fpeAmgSolver>
        addfpeAmgSolverAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
Foam::fpeAmgSolver::fpeAmgSolver
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
    ),
    nSmoothingSteps_(readLabel(dict.lookup("nSmoothingSteps"))),
    mFactor_(readScalar(dict.lookup("mFactor")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::lduSolverPerformance Foam::fpeAmgSolver::solve
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    // Prepare solver performance
    lduSolverPerformance solverPerf(typeName, fieldName());

    scalar normFactor = this->normFactor(x, b, cmpt);

    // Calculate initial residual
    solverPerf.initialResidual() =
        gSumMag(amg_.residual(x, b, cmpt))/normFactor;

    solverPerf.finalResidual() = solverPerf.initialResidual();

    if (!stop(solverPerf))
    {
        scalarField psiSave(x.size());

        do
        {
            amg_.cycle(x, b, cmpt);

            if (solverPerf.nIterations() % nSmoothingSteps_ == 0)
            {
                psiSave = x;
            }
            else if (solverPerf.nIterations() % nSmoothingSteps_ == 1)
            {
                forAll (x, i)
                {
                    x[i] += (mFactor_ - 1)*(x[i] - psiSave[i]);
                }
            }

            // Re-calculate residual
            solverPerf.finalResidual() =
                gSumMag(amg_.residual(x, b, cmpt))/normFactor;
            solverPerf.nIterations()++;
        } while (!stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
