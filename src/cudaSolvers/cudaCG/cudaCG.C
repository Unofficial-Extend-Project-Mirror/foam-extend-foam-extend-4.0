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

\*---------------------------------------------------------------------------*/

#include "cudaCG.H"

// Preconditioner and solver are hardwired due to
// code structure of Cusp library. Below are
// external functions with solver/preconditioner combinations

// CG with diagonal preconditioning
extern "C" void cgDiag
(
    cuspEquationSystem* ces,
    cudaSolverPerformance* solverParam
);

// CG with Ainv preconditioning
extern "C" void cgAinv
(
    cuspEquationSystem* ces,
    cudaSolverPerformance* solverPerf,
    ValueType drop_tolerance,
    int nonzero_per_row,
    bool lin_dropping,
    int lin_param
);

// CG with amg preconditioning
extern "C" void cgAmg
(
    cuspEquationSystem* ces,
    cudaSolverPerformance* solverParam,
    ValueType theta
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cudaCG, 0);

    lduSolver::addsymMatrixConstructorToTable<cudaCG>
        addcudaCGSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// - Construct from matrix and solver data stream
Foam::cudaCG::cudaCG
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

Foam::lduSolverPerformance Foam::cudaCG::solve
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
    cuspEquationSystem ces = createSymCuspMatrix(matrix(), x, b);

    // Call solver externally
    word preconName(dict().lookup("preconditioner"));
    if (preconName == "diagonal")
    {
        cgDiag(&ces, &solverPerf);
    }
    else if(preconName == "Ainv")
    {
        cgAinv
        (
            &ces,
            &solverPerf,
            dict().lookupOrDefault<scalar>("dropTolerance", 0.1),
            dict().lookupOrDefault<label>("nonzeroPerRow", -1),
            dict().lookupOrDefault<bool>("LinDropping", false),
            dict().lookupOrDefault<label>("LinParameter", 1)
        );
    }
    else if(preconName == "amg")
    {
        cgAmg
        (
            &ces,
            &solverPerf,
            dict().lookupOrDefault<scalar>("theta", 0)
        );
    }
    else
    {
        FatalErrorIn("cudaCG::solver()")
            << "Unknown preconditioner name. "
            << "Options are:" << nl
            << "(" << nl
            << "diagonal" << nl
            << "Ainv" << nl
            << "amg" << nl
            << ")" << nl
            << abort(FatalError);
    }

    // Copy the x vector back to Openfoam
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
