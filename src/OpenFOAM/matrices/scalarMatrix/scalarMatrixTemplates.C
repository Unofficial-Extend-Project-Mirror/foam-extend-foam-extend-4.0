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

#include "scalarMatrix.H"
#include "Swap.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::scalarMatrix::solve
(
    Matrix<scalar>& tmpMatrix,
    Field<T>& sourceSol
)
{
    label n = tmpMatrix.n();

    // Elimination
    for (register label i=0; i<n; i++)
    {
        label iMax = i;
        scalar largestCoeff = mag(tmpMatrix[iMax][i]);

        // Swap entries around to find a good pivot
        for (register label j=i+1; j<n; j++)
        {
            if (mag(tmpMatrix[j][i]) > largestCoeff)
            {
                iMax = j;
                largestCoeff = mag(tmpMatrix[iMax][i]);
            }
        }

        if (i != iMax)
        {
            //Info<< "Pivoted on " << i << " " << iMax << endl;

            for (register label k=i; k<n; k++)
            {
                Swap(tmpMatrix[i][k], tmpMatrix[iMax][k]);
            }
            Swap(sourceSol[i], sourceSol[iMax]);
        }

        // Check that the system of equations isn't singular
        if (mag(tmpMatrix[i][i]) < 1e-20)
        {
            FatalErrorIn("scalarMatrix::solve()")
                << "Singular Matrix"
                << exit(FatalError);
        }

        // Reduce to upper triangular form
        for (register label j=i+1; j<n; j++)
        {
            sourceSol[j] -= sourceSol[i]*(tmpMatrix[j][i]/tmpMatrix[i][i]);

            for (register label k=n-1; k>=i; k--)
            {
                tmpMatrix[j][k] -=
                    tmpMatrix[i][k]*tmpMatrix[j][i]/tmpMatrix[i][i];
            }
        }
    }

    // Back-substitution
    for (register label j=n-1; j>=0; j--)
    {
        T ntempvec = pTraits<T>::zero;

        for (register label k=j+1; k<n; k++)
        {
            ntempvec += tmpMatrix[j][k]*sourceSol[k];
        }

        sourceSol[j] = (sourceSol[j] - ntempvec)/tmpMatrix[j][j];
    }
}


template<class T>
void Foam::scalarMatrix::solve(Field<T>& psi, const Field<T>& source) const
{
    Matrix<scalar> tmpMatrix = *this;
    psi = source;
    solve(tmpMatrix, psi);
}


template<class T>
void Foam::scalarMatrix::LUBacksubstitute
(
    const Matrix<scalar>& luMatrix,
    const labelList& pivotIndices,
    Field<T>& sourceSol
)
{
    label n = luMatrix.n();

    label ii = 0;

    for (register label i=0; i<n; i++)
    {
        label ip = pivotIndices[i];
        T sum = sourceSol[ip];
        sourceSol[ip] = sourceSol[i];
        const scalar* __restrict__ luMatrixi = luMatrix[i];

        if (ii != 0)
        {
            for (label j=ii-1; j<i; j++)
            {
                sum -= luMatrixi[j]*sourceSol[j];
            }
        }
        else if (sum != pTraits<T>::zero)
        {
            ii = i+1;
        }

        sourceSol[i] = sum;
    }

    for (register label i=n-1; i>=0; i--)
    {
        T sum = sourceSol[i];
        const scalar* __restrict__ luMatrixi = luMatrix[i];

        for (register label j=i+1; j<n; j++)
        { 
            sum -= luMatrixi[j]*sourceSol[j];
        }

        sourceSol[i] = sum/luMatrixi[i];
    }
}


template<class T>
void Foam::scalarMatrix::LUsolve
(
    Matrix<scalar>& matrix,
    Field<T>& sourceSol
)
{
    labelList pivotIndices(matrix.n());
    LUDecompose(matrix, pivotIndices);
    LUBacksubstitute(matrix, pivotIndices, sourceSol);
}


// ************************************************************************* //
