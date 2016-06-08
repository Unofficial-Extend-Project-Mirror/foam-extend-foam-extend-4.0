// Eigen-decomposition for symmetric 3x3 real matrices.
// Public domain, copied from the public domain Java library JAMA.

#ifndef eig3_H
#define eig3_H

#include "tensor.H"
#include "symmTensor.H"
#include "diagTensor.H"

//#ifndef _eig_h

// Symmetric matrix A => eigenvectors in columns of V,
// corresponding eigenvalues in d.

void eigen_decomposition
(
    const Foam::tensor& A,
    Foam::tensor& V,
    Foam::vector& d
);

void eigen_decomposition
(
    const Foam::symmTensor& sA,
    Foam::tensor& V,
    Foam::vector& d
);

void eigen_decomposition
(
    const Foam::symmTensor& sA,
    Foam::tensor& V,
    Foam::diagTensor& d
);

#endif
