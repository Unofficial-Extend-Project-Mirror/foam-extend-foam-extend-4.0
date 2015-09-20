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

Class
    Foam::cudaSolver

Description
    Base class for GPU based solvers

SourceFiles
    cudaSolver.C

Author
    Dominik Christ, Wikki Ltd.
    Based on Cufflink library by Daniel P. Combest

\*---------------------------------------------------------------------------*/
#ifndef cudaTypes_H
#define cudaTypes_H

#if defined(WM_SP)

#   define MPI_SCALAR MPI_FLOAT
    typedef float ValueType;

#elif defined(WM_DP)

#   define MPI_SCALAR MPI_DOUBLE
    typedef double ValueType;

#elif defined(WM_LDP)

#   error Long double not supported in CUDA

#endif

#include <cusp/coo_matrix.h>
#include <cusp/memory.h>
#include <thrust/host_vector.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef cusp::device_memory MemorySpace;
typedef cusp::host_memory hostMemorySpace;
typedef int IndexType;

// used to prevent floating point exception
const ValueType SMALL = 1e-20;


struct cuspEquationSystem
{
    cusp::coo_matrix<IndexType,ValueType,hostMemorySpace> A;
    cusp::array1d< ValueType, hostMemorySpace> X;
    cusp::array1d< ValueType, hostMemorySpace> B;
    int nCells;
    int nFaces;
};


struct cudaSolverPerformance
{
    int minIter;
    int maxIter;
    double relTol;
    double tol;

    int nIterations;
    double iRes;
    double fRes;
    bool converged;
    bool singular;

    bool debugCusp;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
