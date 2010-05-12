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

#include "GAMGSolver.H"
#include "vector2D.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::GAMGSolver::scalingFactor
(
    scalarField& x,
    const scalarField& b,
    const scalarField& Ax,
    const scalarField& D
) const
{
    scalar scalingFactorNum = 0.0;
    scalar scalingFactorDenom = 0.0;

    forAll(x, i)
    {
        scalingFactorNum += b[i]*x[i];
        scalingFactorDenom += Ax[i]*x[i];

        // While the matrix-multiply done for the scaling it is
        // possible to perform a point-Jacobi smoothing operation cheaply
        x[i] += (b[i] - Ax[i])/D[i];
    }

    vector2D scalingVector(scalingFactorNum, scalingFactorDenom);
    reduce(scalingVector, sumOp<vector2D>());
    return scalingVector.x()/stabilise(scalingVector.y(), VSMALL);
}


Foam::scalar Foam::GAMGSolver::scalingFactor
(
    scalarField& Ax,
    const lduMatrix& A,
    scalarField& x,
    const FieldField<Field, scalar>& coupleLevelBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaceLevel,
    const scalarField& b,
    const direction cmpt
) const
{
    A.Amul
    (
        Ax,
        x,
        coupleLevelBouCoeffs,
        interfaceLevel,
        cmpt
    );

    return scalingFactor
    (
        x,
        b,
        Ax,
        A.diag()
    );
}


// ************************************************************************* //
