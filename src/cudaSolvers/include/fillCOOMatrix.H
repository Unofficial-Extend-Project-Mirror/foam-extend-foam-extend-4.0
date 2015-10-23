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

Please note that cufflink is not approved or endorsed by ESI-OpenCFD®
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

// Fill in the rest of the diag (rows and col),
// upper and upper.rows and upper.cols
cusp::coo_matrix<IndexType, ValueType, MemorySpace> A(ces->A);
cusp::array1d<ValueType, MemorySpace> X(ces->X);
cusp::array1d<ValueType, MemorySpace> B(ces->B);

// Determine row indices of diagonal values and fill A COO matrix
thrust::sequence(A.row_indices.begin(), A.row_indices.begin() + ces->nCells);

// Determine column indices of diagonal values and fill A COO matrix
thrust::sequence
(
    A.column_indices.begin(),
    A.column_indices.begin() + ces->nCells
);

// Copy values of lower into upper in COO matrix
thrust::copy
(
    A.values.begin() + ces->nCells,
    A.values.begin() + ces->nCells + ces->nFaces,
    A.values.begin() + ces->nCells + ces->nFaces
);

// Copy row indices of lower to columns of upper into A COO matrix
thrust::copy
(
    A.row_indices.begin() + ces->nCells,
    A.row_indices.begin() + ces->nCells + ces->nFaces,
    A.column_indices.begin() + ces->nCells + ces->nFaces
);

// Copy column indices of lower to rows of upper into A COO matrix
thrust::copy
(
    A.column_indices.begin() + ces->nCells,
    A.column_indices.begin() + ces->nCells + ces->nFaces,
    A.row_indices.begin() + ces->nCells + ces->nFaces
);

A.sort_by_row_and_column(); // Speeds code up a little bit more
