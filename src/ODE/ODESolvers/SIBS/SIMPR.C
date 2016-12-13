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

#include "SIBS.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::SIBS::SIMPR
(
    const scalar xStart,
    const scalarField& y,
    const scalarField& dydx,
    const scalarField& dfdx,
    const scalarSquareMatrix& dfdy,
    const scalar deltaX,
    const label nSteps,
    scalarField& yEnd
) const
{
    scalar h = deltaX/nSteps;
    const label nEqns = ode_.nEqns();

    scalarSquareMatrix a(nEqns);
    for (register label i=0; i<nEqns; i++)
    {
        for (register label j=0; j<nEqns; j++)
        {
            a[i][j] = -h*dfdy[i][j];
        }
        ++a[i][i];
    }

    labelList pivotIndices(nEqns);
    scalarSquareMatrix::LUDecompose(a, pivotIndices);

    for (register label i=0; i<nEqns; i++)
    {
        yEnd[i] = h*(dydx[i] + h*dfdx[i]);
    }

    scalarSquareMatrix::LUBacksubstitute(a, pivotIndices, yEnd);

    scalarField del(yEnd);
    scalarField ytemp(nEqns);

    for (register label i=0; i<nEqns; i++)
    {
        ytemp[i] = y[i] + del[i];
    }

    scalar x = xStart + h;

    ode_.derivatives(x, ytemp, yEnd);

    for (register label nn=2; nn<=nSteps; nn++)
    {
        for (register label i=0; i<nEqns; i++)
        {
            yEnd[i] = h*yEnd[i] - del[i];
        }

        simpleMatrix<scalar>::LUBacksubstitute(a, pivotIndices, yEnd);

        for (register label i=0; i<nEqns; i++)
        {
            ytemp[i] += (del[i] += 2.0*yEnd[i]);
        }

        x += h;

        ode_.derivatives(x, ytemp, yEnd);
    }
    for (register label i=0; i<nEqns; i++)
    {
        yEnd[i] = h*yEnd[i] - del[i];
    }

    simpleMatrix<scalar>::LUBacksubstitute(a, pivotIndices, yEnd);

    for (register label i=0; i<nEqns; i++)
    {
        yEnd[i] += ytemp[i];
    }
}


// ************************************************************************* //
