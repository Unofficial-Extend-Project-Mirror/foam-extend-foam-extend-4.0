/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "scalarMatrices.H"
#include "SVD.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::multiply
(
    scalarRectangularMatrix& ans,         // value changed in return
    const scalarRectangularMatrix& A,
    const scalarRectangularMatrix& B
)
{
    if (A.m() != B.n())
    {
        FatalErrorIn
        (
            "multiply("
            "scalarRectangularMatrix& answer "
            "const scalarRectangularMatrix& A, "
            "const scalarRectangularMatrix& B)"
        )   << "A and B must have identical inner dimensions but A.m = "
            << A.m() << " and B.n = " << B.n()
            << abort(FatalError);
    }

    ans = scalarRectangularMatrix(A.n(), B.m(), scalar(0));

    for(register label i = 0; i < A.n(); i++)
    {
        for(register label j = 0; j < B.m(); j++)
        {
            for(register label l = 0; l < B.n(); l++)
            {
                ans[i][j] += A[i][l]*B[l][j];
            }
        }
    }
}


void Foam::multiply
(
    scalarRectangularMatrix& ans,         // value changed in return
    const scalarRectangularMatrix& A,
    const scalarRectangularMatrix& B,
    const scalarRectangularMatrix& C
)
{
    if (A.m() != B.n())
    {
        FatalErrorIn
        (
            "multiply("
            "const scalarRectangularMatrix& A, "
            "const scalarRectangularMatrix& B, "
            "const scalarRectangularMatrix& C, "
            "scalarRectangularMatrix& answer)"
        )   << "A and B must have identical inner dimensions but A.m = "
            << A.m() << " and B.n = " << B.n()
            << abort(FatalError);
    }

    if (B.m() != C.n())
    {
        FatalErrorIn
        (
            "multiply("
            "const scalarRectangularMatrix& A, "
            "const scalarRectangularMatrix& B, "
            "const scalarRectangularMatrix& C, "
            "scalarRectangularMatrix& answer)"
        )   << "B and C must have identical inner dimensions but B.m = "
            << B.m() << " and C.n = " << C.n()
            << abort(FatalError);
    }

    ans = scalarRectangularMatrix(A.n(), C.m(), scalar(0));

    for(register label i = 0; i < A.n(); i++)
    {
        for(register label g = 0; g < C.m(); g++)
        {
            for(register label l = 0; l < C.n(); l++)
            {
                scalar ab = 0;
                for(register label j = 0; j < A.m(); j++)
                {
                    ab += A[i][j]*B[j][l];
                }
                ans[i][g] += C[l][g] * ab;
            }
        }
    }
}


void Foam::multiply
(
    scalarRectangularMatrix& ans,         // value changed in return
    const scalarRectangularMatrix& A,
    const DiagonalMatrix<scalar>& B,
    const scalarRectangularMatrix& C
)
{
    if (A.m() != B.size())
    {
        FatalErrorIn
        (
            "multiply("
            "const scalarRectangularMatrix& A, "
            "const DiagonalMatrix<scalar>& B, "
            "const scalarRectangularMatrix& C, "
            "scalarRectangularMatrix& answer)"
        )   << "A and B must have identical inner dimensions but A.m = "
            << A.m() << " and B.n = " << B.size()
            << abort(FatalError);
    }

    if (B.size() != C.n())
    {
        FatalErrorIn
        (
            "multiply("
            "const scalarRectangularMatrix& A, "
            "const DiagonalMatrix<scalar>& B, "
            "const scalarRectangularMatrix& C, "
            "scalarRectangularMatrix& answer)"
        )   << "B and C must have identical inner dimensions but B.m = "
            << B.size() << " and C.n = " << C.n()
            << abort(FatalError);
    }

    ans = scalarRectangularMatrix(A.n(), C.m(), scalar(0));

    for(register label i = 0; i < A.n(); i++)
    {
        for(register label g = 0; g < C.m(); g++)
        {
            for(register label l = 0; l < C.n(); l++)
            {
                ans[i][g] += C[l][g] * A[i][l]*B[l];
            }
        }
    }
}


Foam::RectangularMatrix<Foam::scalar> Foam::SVDinv
(
    const scalarRectangularMatrix& A,
    scalar minCondition
)
{
    SVD svd(A, minCondition);
    return svd.VSinvUt();
}


// ************************************************************************* //
