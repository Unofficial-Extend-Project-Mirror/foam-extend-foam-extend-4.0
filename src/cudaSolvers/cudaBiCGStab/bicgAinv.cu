/**********************************************************************\
  ______  __    __   _______  _______  __       __   __   __   __  ___
 /      ||  |  |  | |   ____||   ____||  |     |  | |  \ |  | |  |/  /
|  ,----'|  |  |  | |  |__   |  |__   |  |     |  | |   \|  | |  '  /
|  |     |  |  |  | |   __|  |   __|  |  |     |  | |  . `  | |    <
|  `----.|  `--'  | |  |     |  |     |  `----.|  | |  |\   | |  .  \
 \______| \______/  |__|     |__|     |_______||__| |__| \__| |__|\__\

Cuda For FOAM Link

cufflink is a library for linking numerical methods based on Nvidia's
Compute Unified Device Architecture (CUDA™) C/C++ programming language
and OpenFOAM®.

Please note that cufflink is not approved or endorsed by OpenCFD®
Limited, the owner of the OpenFOAM® and OpenCFD® trademarks and
producer of OpenFOAM® software.

The official web-site of OpenCFD® Limited is www.openfoam.com .

------------------------------------------------------------------------
This file is part of cufflink.

    cufflink is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    cufflink is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cufflink.  If not, see <http://www.gnu.org/licenses/>.

    Author
    Daniel P. Combest.  All rights reserved.
    Modifications by Dominik Christ, Wikki Ltd.

    Description
    diagonal preconditioned conjugate gradient
    solver for symmetric Matrices using a CUSP CUDA™ based solver.

\**********************************************************************/

#include "cudaTypes.H"

// CUSP Includes
#include <cusp/detail/config.h>
#include <cusp/verify.h>
#include <cusp/precond/ainv.h>
#include <cusp/coo_matrix.h>
#include <cusp/blas.h>
#include <cusp/multiply.h>

extern "C" void bicgAinv
(
    cuspEquationSystem* ces,
    cudaSolverPerformance* solverPerf,
    ValueType drop_tolerance,
    int nonzero_per_row,
    bool lin_dropping,
    int lin_param
)
{
    cusp::coo_matrix<IndexType, ValueType, MemorySpace> A(ces->A);
    cusp::array1d<ValueType, MemorySpace> X(ces->X);
    cusp::array1d<ValueType, MemorySpace> B(ces->B);

    // Fill in the rest of the diag (rows and col)
    // Determine row indices of diagonal values and fill A COO matrix
    thrust::sequence
    (
        A.row_indices.begin(),
        A.row_indices.begin() + ces->nCells
    );

    // Determine column indices of diagonal values and fill A COO matrix
    thrust::sequence
    (
        A.column_indices.begin(),
        A.column_indices.begin() + ces->nCells
    );

    // Sorted coo by row and column. speeds code up a little bit more
    A.sort_by_row_and_column();

    // Set storage
    // #include "cufflink/setGPUStorage.H"

    if (solverPerf->debugCusp)
    {
        std::cout << "   Using Ainv preconditioner\n";
    }

    cusp::precond::bridson_ainv<ValueType, MemorySpace> M
    (
         A,
         drop_tolerance,
         nonzero_per_row,
         lin_dropping,
         lin_param
    );

    // Start Krylov solver
    assert(A.num_rows == A.num_cols);   // Sanity check

    const size_t N = A.num_rows;

    // Allocate workspace
    cusp::array1d<ValueType,MemorySpace> y(N);

    cusp::array1d<ValueType,MemorySpace> p(N);
    cusp::array1d<ValueType,MemorySpace> r(N);
    cusp::array1d<ValueType,MemorySpace> r_star(N);
    cusp::array1d<ValueType,MemorySpace> s(N);
    cusp::array1d<ValueType,MemorySpace> Mp(N);
    cusp::array1d<ValueType,MemorySpace> AMp(N);
    cusp::array1d<ValueType,MemorySpace> Ms(N);
    cusp::array1d<ValueType,MemorySpace> AMs(N);

    // y <- Ax
    cusp::multiply(A, X, y);

    // Define the normalization factor
    ValueType normFactor = 1.0;

#   include "buildNormFactor.H"

    // r <- b - A*x
    cusp::blas::axpby(B, y, r, ValueType(1), ValueType(-1));

    // p <- r
    cusp::blas::copy(r, p);

    // r_star <- r
    cusp::blas::copy(r, r_star);

    ValueType r_r_star_old = cusp::blas::dotc(r_star, r);


    ValueType normR = cusp::blas::nrm2(r)/normFactor;
    ValueType normR0 = normR;   // Initial residual
    solverPerf->iRes    = normR0;
    int count = 0;

    while
    (
        normR > (solverPerf->tol)
     && count < (solverPerf->maxIter)
     && normR/normR0 >= (solverPerf->relTol)
     || count < solverPerf->minIter
    )
    {
        // Mp = M*p
        cusp::multiply(M, p, Mp);

        // AMp = A*Mp
        cusp::multiply(A, Mp, AMp);

        // alpha = (r_j, r_star) / (A*M*p, r_star)
        ValueType alpha = r_r_star_old/cusp::blas::dotc(r_star, AMp);

        // s_j = r_j - alpha * AMp
        cusp::blas::axpby(r, AMp, s, ValueType(1), ValueType(-alpha));

        ValueType normS = cusp::blas::nrm2(s)/normFactor;

        if
        (
           !(
                normS > (solverPerf->tol)
             && normS/normR0 >= (solverPerf->relTol)
            )
        )
        {
            // is this right?
            // x += alpha*M*p_j
            cusp::blas::axpby(X, Mp, X, ValueType(1), ValueType(alpha));

            // y <- Ax
            cusp::multiply(A, X, y);

            // r <- b - A*x
            cusp::blas::axpby(B, y, r, ValueType(1), ValueType(-1));
            // not sure we should be measuring r here...need to look again.
            normR = cusp::blas::nrm2(r)/normFactor;

            count++;

            break;
        }

        // Ms = M*s_j
        cusp::multiply(M, s, Ms);

        // AMs = A*Ms
        cusp::multiply(A, Ms, AMs);

        // omega = (AMs, s) / (AMs, AMs)
        ValueType omega = cusp::blas::dotc(AMs, s)/cusp::blas::dotc(AMs, AMs);

        // x_{j+1} = x_j + alpha*M*p_j + omega*M*s_j
        cusp::blas::axpbypcz(X, Mp, Ms, X, ValueType(1), alpha, omega);

        // r_{j+1} = s_j - omega*A*M*s
        cusp::blas::axpby(s, AMs, r, ValueType(1), -omega);

        // beta_j = (r_{j+1}, r_star) / (r_j, r_star) * (alpha/omega)
        ValueType r_r_star_new = cusp::blas::dotc(r_star, r);
        ValueType beta = (r_r_star_new/r_r_star_old)*(alpha/omega);
        r_r_star_old = r_r_star_new;

        // p_{j+1} = r_{j+1} + beta*(p_j - omega*A*M*p)
        cusp::blas::axpbypcz(r, p, AMp, p, ValueType(1), beta, -beta*omega);

        // not dure we should be measuring r here...need to look again.
        normR = cusp::blas::nrm2(r)/normFactor;

        count++;
    }
    // End Krylov Solver

    // Final residual
    solverPerf->fRes = cusp::blas::nrm2(r)/normFactor;
    solverPerf->nIterations = count;

    // converged?
    if
    (
        solverPerf->fRes<=solverPerf->tol
     || solverPerf->fRes/solverPerf->iRes <= solverPerf->relTol
    )
    {
        solverPerf->converged = true;
    }
    else
    {
        solverPerf->converged = false;
    }

    // Pass the solution vector back
    ces->X = X;
}
