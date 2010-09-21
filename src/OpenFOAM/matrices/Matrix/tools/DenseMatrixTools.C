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

Class
    DenseMatrixTools

\*---------------------------------------------------------------------------*/

#include "DenseMatrixTools.H"
#include "Swap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Solve with pivoting.  Note: Destroys matrix and b
template<class T>
static void Foam::DenseMatrixTools::solve
(
    SquareMatrix<T>& A,
    List<T>& x,
    List<T>& b
)
{
    const label nRows = A.n();

    // Create and initialise permutation array
    labelList p(nRows);

    forAll (p, i)
    {
        p[i] = i;
    }


    // Eliminate equations
    for (label k = 0; k < nRows - 1; k++)
    {
        // Swap rows
        label m = k;

        for(label i = k + 1; i < nRows; ++i)
        {
            if (Foam::mag(A[p[i]][k]) > Foam::mag(A[p[m]][k]))
            {
                m = i;
            }

            Swap(p[k], p[m]);
        }


        for (label i = k + 1; i < nRows; i++)
        {
            label pi = p[i];
            label pk = p[k];

//             T r = A[pi][k]/A[pk][k];
            T r = cmptDivide(A[pi][k], A[pk][k]);

            for (label j = k + 1; j < nRows; j++)
            {
//                 A[pi][j] -= A[pk][j]*r;
                A[pi][j] -= cmptMultiply(A[pk][j], r);
            }

//             b[pi] -= b[pk]*r;
            b[pi] -= cmptMultiply(b[pk], r);
        }
    }

    // Back substitute
    for (label i = nRows - 1; i >= 0; i--)
    {
        label pi = p[i];

        T sum = b[pi];

        for (label j = i + 1; j < nRows; j++)
        {
//             sum -= A[pi][j]*x[j];
            sum -= cmptMultiply(A[pi][j], x[j]);
        }

//         x[i] = sum/A[pi][i];
        x[i] = cmptDivide(sum, A[pi][i]);
    }
}


template<class Form, class Type>
static void Foam::DenseMatrixTools::qrDecompose
(
    const label nCols,
    FieldField<Field, Type>& A,
    Matrix<Form, Type>& R
)
{
    // Note: consider Arnoldi algorithm for speed-up.  HJ, 14/Sep/2006

    for (label j = 0; j < nCols; j++)
    {
        R[j][j] = Foam::sqrt(gSumSqr(A[j]));

        if (mag(R[j][j]) > SMALL)
        {
            A[j] /= R[j][j];
        }

        for (label k = j + 1; k < nCols; k++)
        {
            const Field<Type>& Aj = A[j];
            Field<Type>& Ak = A[k];
            Type& Rjk = R[j][k];

            Rjk = gSumProd(Aj, Ak);

            for (label i = 0; i < Ak.size(); i++)
            {
                Ak[i] -= Rjk*Aj[i];
            }
        }
    }
}


// ************************************************************************* //
