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

Description
    Preconditioned Conjugate Gradient solver with run-time selectable
    preconditioning

Author
    Dominik Christ, Wikki Ltd.
    Based on Cufflink library by Daniel P. Combest

\*---------------------------------------------------------------------------*/

#include "cudaBiCGStab.H"

// Preconditioner and solver are hardwired due to
// code structure of Cusp library. Below are
// external functions with solver/preconditioner combinations

// BiCG with diagonal preconditioning
extern "C" void bicgDiag
(
    cuspEquationSystem* ces,
    cudaSolverPerformance* solverParam
);

// BiCG with Ainv preconditioning
extern "C" void bicgAinv
(
    cuspEquationSystem* ces,
    cudaSolverPerformance* solverPerf,
    ValueType drop_tolerance,
    int nonzero_per_row,
    bool lin_dropping,
    int lin_param
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cudaBiCGStab, 0);

    lduSolver::addasymMatrixConstructorToTable<cudaBiCGStab>
        addcudaBiCGStabAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// - Construct from matrix and solver data stream
Foam::cudaBiCGStab::cudaBiCGStab
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& dict
)
:
    cudaSolver
    (
        fieldName,
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces,
        dict
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::lduSolverPerformance Foam::cudaBiCGStab::solve
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    // Initialize Cusp solver perfomance
    cudaSolverPerformance solverPerf = cudaSolverPerformanceDefault();
    solverPerf.minIter = minIter(); // Minimum iterations
    solverPerf.maxIter = maxIter(); // Maximum iterations
    solverPerf.relTol = relTolerance(); // Relative tolerance
    solverPerf.tol = tolerance(); // Tolerance

    if (lduMatrix::debug >= 2)
    {
         solverPerf.debugCusp = true;
    }

    // Initialize and copy matrix data to GPU
    cuspEquationSystem ces = createAsymCuspMatrix(matrix(), x, b);

    // Call solver externally
    word preconName(dict().lookup("preconditioner"));
    if (preconName == "diagonal")
    {
        bicgDiag(&ces, &solverPerf);
    }
    else if(preconName == "Ainv")
    {
        bicgAinv
        (
            &ces,
            &solverPerf,
            dict().lookupOrDefault<scalar>("dropTolerance", 0.1),
            dict().lookupOrDefault<label>("nonzeroPerRow", -1),
            dict().lookupOrDefault<bool>("LinDropping", false),
            dict().lookupOrDefault<label>("LinParameter", 1)
        );
    }
    else
    {
        FatalErrorIn("cudaBiCGStab::solver()")
            << "Unknown preconditioner name. "
            << "Options are:" << nl
            << "(" << nl
            << "diagonal" << nl
            << "Ainv" << nl
            << ")" << nl
            << abort(FatalError);
    }

    // copy the x vector back to Openfoam
    thrust::copy(ces.X.begin(), ces.X.end(), x.begin());

    // Return solver output
    return lduSolverPerformance
    (
        typeName,
        fieldName(),
        solverPerf.iRes,
        solverPerf.fRes,
        solverPerf.nIterations,
        solverPerf.converged,
        solverPerf.singular
    );

}


// ************************************************************************* //
