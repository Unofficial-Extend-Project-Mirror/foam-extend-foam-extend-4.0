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

#include "splineInterpolateXY.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char*
Foam::NamedEnum<Foam::splineInterpolateXY::splineBCType, 2>::names[] =
{
    "not-a-knot",
    "natural"
};

const Foam::NamedEnum<Foam::splineInterpolateXY::splineBCType, 2>
    Foam::splineInterpolateXY::splineBCTypeNames;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::splineInterpolateXY::splineInterpolateXY
(
    const scalarField& x,
    const scalarField& y,
    const word& startBC,
    const word& endBC
)
:
    X_(x),
    Y_(y),
    Y2_()
{
    setData(splineBCTypeNames[startBC], splineBCTypeNames[endBC]);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::splineInterpolateXY::setData
(
    const splineBCType& startBC,
    const splineBCType& endBC
)
{
    // Objective now is to approximate second derivative at given data points

    const scalar n = X_.size();

    // Define matrix
    simpleMatrix<scalar> A(n);

    // Fill matrix element by element
    for (label j = 1; j < n - 1; j++)
    {
        A.source()[j] = (Y_[j+1] - Y_[j])/(X_[j + 1] - X_[j])
            - (Y_[j] - Y_[j - 1])/(X_[j] - X_[j - 1]);

        A[j][j - 1] = (X_[j] - X_[j - 1])/6;
        A[j][j] = (X_[j + 1] - X_[j - 1])/3;
        A[j][j + 1] = (X_[j + 1] - X_[j])/6;
    }

    // Explicitly set boundary values by adjusting
    // the first and last row of the matrix
    // (which also make the matrix a square)
    scalarField firstRow(n, 0);

    // Start BC
    if (startBC == NOT_A_KNOT)
    {
        // Not-a-knot takes the same third derivative
        // as in the second and second last data point
        // matrix-wise this is:
        firstRow[0] = -1./6;
        firstRow[1] = 1./3;
        firstRow[2] = -1./6;
    }
    else if (startBC == NATURAL)
    {
        // Natural spline assumes zero second derivative
        // in end knots, i.e. Y2_[0] = 0 and Y2_[n] = 0
        firstRow[0] = 1;
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::splineInterpolateXY::setData\n"
            "(\n"
            "    const splineBCType& startBC,\n"
            "    const splineBCType& endBC\n"
            ")"
        )   << "Unknown splineBC: " << splineBCTypeNames[startBC]
            << ". Should be " << splineBCTypeNames
            << abort(FatalError);
    }

    scalarField lastRow(n, 0);

    // End BC
    if (endBC == NOT_A_KNOT)
    {
        // Not-a-knot takes the same third derivative
        // as in the second and second last data point
        // matrix-wise this is:

        lastRow[n - 3] = -1./6;
        lastRow[n - 2] = 1./3;
        lastRow[n - 1] = -1./6;

    }
    else if (endBC == NATURAL)
    {
        // Natural spline assumes zero second derivative
        // in end knots, i.e. Y2_[0] = 0 and Y2_[n] = 0
        lastRow[n - 1] = 1;
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::splineInterpolateXY::setData\n"
            "(\n"
            "    const splineBCType& startBC,\n"
            "    const splineBCType& endBC\n"
            ")"
        )   << "Unknown splineBC: " << splineBCTypeNames[endBC]
            << ". Should be " << splineBCTypeNames
            << abort(FatalError);
    }

    // Set the first row
    forAll (firstRow, i)
    {
        A[0][i] = firstRow[i];
    }
    A.source()[0] = 0;

    // Set the last row
    forAll(lastRow, i)
    {
        A[n-1][i] = lastRow[i];
    }
    A.source()[n-1] = 0;

    // Solve system and obtain second derivatives
    Y2_ = A.LUsolve();
}


Foam::scalar Foam::splineInterpolateXY::interpolate
(
    const scalar XStar
) const
{
    // Find indices of data points that enclose XStar using bisection
    label nLow = 0;
    label nHigh = X_.size() - 1;

    while ((nHigh - nLow) > 1)
    {
        label n = (nHigh + nLow)/2;
        if (X_[n] > XStar)
        {
            nHigh = n;
        }
        else
        {
            nLow = n;
        }
    }

    // Set up coefficients
    scalar h = X_[nHigh] - X_[nLow];
    scalar a = (X_[nHigh] - XStar)/(h + SMALL);
    scalar b = (XStar - X_[nLow])/(h + SMALL);

    // Svaluate spline
    return a*Y_[nLow] + b*Y_[nHigh]
        + ((pow3(a) - a)*Y2_[nLow] + (pow3(b) - b)*Y2_[nHigh])*sqr(h)/6;
}


// ************************************************************************* //
