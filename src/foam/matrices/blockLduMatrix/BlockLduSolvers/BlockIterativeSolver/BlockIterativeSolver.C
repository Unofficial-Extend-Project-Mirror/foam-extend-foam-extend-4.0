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

Description
    Virtual base class for block iterative solvers

\*---------------------------------------------------------------------------*/

#include "BlockIterativeSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
template<class Type>
Foam::BlockIterativeSolver<Type>::BlockIterativeSolver
(
    const word& fieldName,
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
:
    BlockLduSolver<Type>(fieldName, matrix, dict),
    tolerance_(readScalar(this->dict().lookup("tolerance"))),
    relTolerance_(readScalar(this->dict().lookup("relTol"))),
    minIter_(readLabel(this->dict().lookup("minIter"))),
    maxIter_(readLabel(this->dict().lookup("maxIter")))
{}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::BlockIterativeSolver<Type>::normFactor
(
    Field<Type>& x,
    const Field<Type>& b
) const
{
    const BlockLduMatrix<Type>& matrix = this->matrix_;

    // Calculate the normalisation factor
    const label nRows = x.size();

    Field<Type> pA(nRows);
    Field<Type> wA(nRows);

    // Calculate reference value of x
    Type xRef = gAverage(x);

    // Calculate A.x
    matrix.Amul(wA, x);

    // Calculate A.xRef, temporarily using pA for storage
    matrix.Amul
    (
        pA,
        Field<Type>(nRows, xRef)
    );

    scalar normFactor = gSum(mag(wA - pA) + mag(b - pA)) + this->small_;

    if (BlockLduMatrix<Type>::debug >= 2)
    {
        Info<< "Iterative solver normalisation factor = "
            << normFactor << endl;
    }

    return normFactor;
}


template<class Type>
bool Foam::BlockIterativeSolver<Type>::stop
(
    BlockSolverPerformance<Type>& solverPerf
) const
{
    if (solverPerf.nIterations() < minIter_)
    {
        return false;
    }

    if
    (
        solverPerf.nIterations() >= maxIter_
     || solverPerf.checkConvergence(tolerance_, relTolerance_)
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
