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

\*---------------------------------------------------------------------------*/

#include "cudaSolver.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cudaSolver, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cudaSolver::cudaSolver
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
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

cuspEquationSystem Foam::cudaSolver::createSymCuspMatrix
(
    const lduMatrix& matrix,
    const scalarField& x,
    const scalarField& b
) const
{
    cuspEquationSystem ces;
    ces.nCells = x.size();
    ces.nFaces = matrix.lower().size();

    ces.A = cusp::coo_matrix<IndexType, ValueType, hostMemorySpace>
    (
        ces.nCells,
        ces.nCells,
        ces.nCells+2*ces.nFaces
    );
    ces.X = cusp::array1d< ValueType, hostMemorySpace>(ces.nCells);
    ces.B = cusp::array1d< ValueType, hostMemorySpace>(ces.nCells);

    // Copy values of lduMatrix diag to A COO matrix
    thrust::copy
    (
        matrix.diag().begin(),
        matrix.diag().end(),
        ces.A.values.begin()
    );

    // Copy values of lduMatrix lower to A COO matrix
    thrust::copy
    (
        matrix.lower().begin(),
        matrix.lower().end(),
        ces.A.values.begin()+ces.nCells
    );

    // Since matrix is symmetric, do not copy upper values to save time
    // -> performed on GPU
    // Copy row indices of lower into A COO matrix
    thrust::copy
    (
        matrix.lduAddr().upperAddr().begin(),
        matrix.lduAddr().upperAddr().end(),
        ces.A.row_indices.begin()+ces.nCells
    );
    // Copy column indices of lower into A COO matrix
    thrust::copy
    (
        matrix.lduAddr().lowerAddr().begin(),
        matrix.lduAddr().lowerAddr().end(),
        ces.A.column_indices.begin()+ces.nCells
    );

    // Do not initialize the row and column values of diag to save time
    // -> performed on GPU
    // Copy x of lower into x vector
    thrust::copy(x.begin(), x.end(), ces.X.begin());
    // Copy b of lower into b vector
    thrust::copy(b.begin(), b.end(), ces.B.begin());

    return ces;
}

cuspEquationSystem Foam::cudaSolver::createAsymCuspMatrix
(
    const lduMatrix& matrix,
    const scalarField& x,
    const scalarField& b
) const
{
    cuspEquationSystem ces;
    ces.nCells = x.size();
    ces.nFaces = matrix.lower().size();

    ces.A = cusp::coo_matrix<IndexType, ValueType, hostMemorySpace>
    (
        ces.nCells,
        ces.nCells,
        ces.nCells + 2*ces.nFaces
    );
    ces.X = cusp::array1d< ValueType, hostMemorySpace>(ces.nCells);
    ces.B = cusp::array1d< ValueType, hostMemorySpace>(ces.nCells);

    // Copy values from the lduMatrix to our equation system
    // Copy values of lduMatrix diag to A COO matrix
    thrust::copy
    (
        matrix.diag().begin(),
        matrix.diag().end(),
        ces.A.values.begin()
    );

    // Copy values of lduMatrix lower to A COO matrix
    thrust::copy
    (
        matrix.lower().begin(),
        matrix.lower().end(),
        ces.A.values.begin() + ces.nCells
    );

    // Copy values of lduMatrix upper to A COO matrix
    thrust::copy
    (
        matrix.upper().begin(),
        matrix.upper().end(),
        ces.A.values.begin() + ces.nCells + ces.nFaces
    );

    // Copy row and column indices of lower and upper to our equations system
    // do not initialize the row and column values of diag to save time
    // -> performed on GPU

    // Copy row indices of lower into A COO matrix
    thrust::copy
    (
        matrix.lduAddr().upperAddr().begin(),
        matrix.lduAddr().upperAddr().end(),
        ces.A.row_indices.begin() + ces.nCells
    );

    // Copy column indices of lower into A COO matrix
    thrust::copy
    (
        matrix.lduAddr().lowerAddr().begin(),
        matrix.lduAddr().lowerAddr().end(),
        ces.A.column_indices.begin() + ces.nCells
    );

    // Copy row indices of upper into A COO matrix
    thrust::copy
    (
        matrix.lduAddr().lowerAddr().begin(),
        matrix.lduAddr().lowerAddr().end(),
        ces.A.row_indices.begin() + ces.nCells + ces.nFaces
    );

    // Copy column indices of upper into A COO matrix
    thrust::copy
    (
        matrix.lduAddr().upperAddr().begin(),
        matrix.lduAddr().upperAddr().end(),
        ces.A.column_indices.begin() + ces.nCells + ces.nFaces
    );

    // Copy x of lower into x vector
    thrust::copy(x.begin(), x.end(), ces.X.begin());
    // Copy b of lower into b vector
    thrust::copy(b.begin(), b.end(), ces.B.begin());
    return ces;
}

cudaSolverPerformance Foam::cudaSolver::cudaSolverPerformanceDefault() const
{
    cudaSolverPerformance csp;

    csp.minIter = 0;
    csp.maxIter = 1000;
    csp.relTol = 0;
    csp.tol = 1e-6;

    csp.nIterations = 0;
    csp.iRes = -1;
    csp.fRes = -1;
    csp.converged = false;
    csp.singular = false;

    csp.debugCusp = false;

    return csp;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
