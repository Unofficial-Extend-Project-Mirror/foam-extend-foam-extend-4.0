/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    fineBlockAmgLevel

Description
    Finest AMG level container, using matrix and field references

Author
    Klas Jareteg, 2013-04-15

\*---------------------------------------------------------------------------*/

#include "fineBlockAmgLevel.H"
#include "coarseBlockAmgLevel.H"
#include "BlockSolverPerformance.H"
#include "BlockCoeffNorm.H"
#include "BlockCoeffTwoNorm.H"
#include "BlockBiCGStabSolver.H"
#include "BlockCGSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Type>
Foam::fineBlockAmgLevel<Type>::fineBlockAmgLevel
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict,
    const word& coarseningType,
    const label groupSize,
    const label minCoarseEqns,
    const word& smootherType
)
:
    matrix_(matrix),
    dict_(dict),
    coarseningPtr_
    (
        BlockMatrixCoarsening<Type>::New
        (
            coarseningType,
            matrix_,
            dict_,
            groupSize,
            minCoarseEqns
        )
    ),
    smootherPtr_
    (
        BlockLduSmoother<Type>::New
        (
            matrix,
            dict
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::Field<Type>& Foam::fineBlockAmgLevel<Type>::x()
{
    FatalErrorIn("Field<Type>& Foam::fineBlockAmgLevel<Type>::x()")
        << "x is not available."
        << abort(FatalError);

    // Dummy return
    return const_cast<Field<Type>&>(Field<Type>::null());
}


template<class Type>
Foam::Field<Type>& Foam::fineBlockAmgLevel<Type>::b()
{
    FatalErrorIn("Field<Type>& Foam::fineBlockAmgLevel<Type>::b()")
        << "b is not available."
        << abort(FatalError);

    // Dummy return
    return const_cast<Field<Type>&>(Field<Type>::null());
}


template<class Type>
void Foam::fineBlockAmgLevel<Type>::residual
(
    const Field<Type>& x,
    const Field<Type>& b,
    Field<Type>& res
) const
{
    matrix_.Amul
    (
        res,
        x
    );

    // residual = b - Ax
    forAll (b, i)
    {
        res[i] = b[i] - res[i];
    }
}


template<class Type>
void Foam::fineBlockAmgLevel<Type>::restrictResidual
(
    const Field<Type>& x,
    const Field<Type>& b,
    Field<Type>& xBuffer,
    Field<Type>& coarseRes,
    bool preSweepsDone
) const
{
    if (preSweepsDone)
    {
        // Calculate residual
        residual(x, b, xBuffer);
    }

    // Here x != 0.  It is assumed that the buffer will contain the residual
    // if no pre-sweeps have been done.  HJ, 4/Sep/2006
    coarseningPtr_->restrictResidual(xBuffer, coarseRes);
}


template<class Type>
void Foam::fineBlockAmgLevel<Type>::prolongateCorrection
(
    Field<Type>& x,
    const Field<Type>& coarseX
) const
{
    coarseningPtr_->prolongateCorrection(x, coarseX);
}


template<class Type>
void Foam::fineBlockAmgLevel<Type>::smooth
(
    Field<Type>& x,
    const Field<Type>& b,
    const label nSweeps
) const
{
    smootherPtr_->smooth(x, b, nSweeps);
}


template<class Type>
void Foam::fineBlockAmgLevel<Type>::solve
(
    Field<Type>& x,
    const Field<Type>& b,
    const scalar tolerance,
    const scalar relTol
) const
{
    Info<< "Fine level solver"<<endl;

    // Create artificial dictionary for finest solution
    dictionary finestDict;
//     finestDict.add("nDirections", "5");
    finestDict.add("minIter", 1);
    finestDict.add("maxIter", 1000);
    finestDict.add("tolerance", tolerance);
    finestDict.add("relTol", relTol);

    if (matrix_.symmetric())
    {
        finestDict.add("preconditioner", "Cholesky");

        BlockSolverPerformance<Type> coarseSolverPerf =
            BlockCGSolver<Type>
            (
                "topLevelCorr",
                matrix_,
                finestDict
            ).solve(x, b);

        if (BlockLduMatrix<Type>::debug >= 2)
        {
            coarseSolverPerf.print();
        }
    }
    else
    {
        finestDict.add("preconditioner", "Cholesky");

        BlockSolverPerformance<Type> coarseSolverPerf =
            BlockBiCGStabSolver<Type>
            (
                "topLevelCorr",
                matrix_,
                finestDict
            ).solve(x, b);

        if (BlockLduMatrix<Type>::debug >= 2)
        {
            coarseSolverPerf.print();
        }
    }
}


template<class Type>
void Foam::fineBlockAmgLevel<Type>::scaleX
(
    Field<Type>& x,
    const Field<Type>& b,
    Field<Type>& xBuffer
) const
{

    // KRJ: 2013-02-05: Removed subfield, creating a new field
    Field<Type> Ax(x.size());

    matrix_.Amul
    (
        reinterpret_cast<Field<Type>&>(Ax),
        x
    );

    scalar scalingFactorNum = sumProd(x, b);
    scalar scalingFactorDenom = sumProd(x, Ax);

    vector scalingVector(scalingFactorNum, scalingFactorDenom, 0);
    reduce(scalingVector, sumOp<vector>());

    // Scale x
    if
    (
        scalingVector[0]*scalingVector[1] <= 0
     || mag(scalingVector[0]) < mag(scalingVector[1])
    )
    {
        // Factor = 1.0, no scaling
    }
    else if (mag(scalingVector[0]) > 2*mag(scalingVector[1]))
    {
        // Max factor = 2
        x *= 2.0;
    }
    else
    {
        // Regular scaling
        x *= scalingVector[0]/stabilise(scalingVector[1], SMALL);
    }
}


template<class Type>
Foam::autoPtr<Foam::BlockAmgLevel<Type> >
Foam::fineBlockAmgLevel<Type>::makeNextLevel() const
{
    if (coarseningPtr_->coarsen())
    {
        return autoPtr<Foam::BlockAmgLevel<Type> >
        (
            new coarseBlockAmgLevel<Type>
            (
                coarseningPtr_->restrictMatrix(),
                dict(),
                coarseningPtr_->type(),
                coarseningPtr_->groupSize(),
                coarseningPtr_->minCoarseEqns(),
                smootherPtr_->type()
            )
        );
    }
    else
    {
        // Final level: cannot coarsen
        return autoPtr<Foam::BlockAmgLevel<Type> >();
    }
}


// ************************************************************************* //
