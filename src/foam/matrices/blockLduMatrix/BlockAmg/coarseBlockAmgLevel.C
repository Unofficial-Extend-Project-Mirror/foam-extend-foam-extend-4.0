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
    coarseBlockAmgLevel

Description
    Coarse AMG level stores matrix, x and b locally, for BlockLduMatrix

Author
    Klas Jareteg, 2012-12-13

\*---------------------------------------------------------------------------*/

#include "coarseBlockAmgLevel.H"
#include "SubField.H"
#include "ICCG.H"
#include "BICCG.H"
#include "vector2D.H"
#include "BlockSolverPerformance.H"
#include "BlockGaussSeidelSolver.H"
#include "BlockBiCGStabSolver.H"
#include "BlockCGSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Type>
Foam::coarseBlockAmgLevel<Type>::coarseBlockAmgLevel
(
    autoPtr<BlockLduMatrix<Type> > matrixPtr,
    const dictionary& dict,
    const word& policyType,
    const label groupSize,
    const label minCoarseEqns,
    const word& smootherType
)
:
    matrixPtr_(matrixPtr),
    x_(matrixPtr_->diag().size(),pTraits<Type>::zero),
    b_(matrixPtr_->diag().size(),pTraits<Type>::zero),
    dict_(dict),
    policyPtr_
    (
        BlockAmgPolicy<Type>::New
        (
            policyType,
            matrixPtr_,
            dict_,
            groupSize,
            minCoarseEqns
        )
    ),
    smootherPtr_
    (
        BlockLduSmoother<Type>::New
        (
            matrixPtr_,
            dict
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::coarseBlockAmgLevel<Type>::~coarseBlockAmgLevel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::Field<Type>& Foam::coarseBlockAmgLevel<Type>::x()
{
    return x_;
}


template<class Type>
Foam::Field<Type>& Foam::coarseBlockAmgLevel<Type>::b()
{
    return b_;
}


template<class Type>
void Foam::coarseBlockAmgLevel<Type>::residual
(
    const Field<Type>& x,
    const Field<Type>& b,
    Field<Type>& res
) const
{
    // Calculate residual
    matrixPtr_->Amul
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
void Foam::coarseBlockAmgLevel<Type>::restrictResidual
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
        // KRJ: 2012-12-14 ::subField removed. Creating buffer locally
        Field<Type> resBuf(x.size());

        Field<Type>& res = reinterpret_cast<Field<Type>&>(resBuf);

        residual(x, b, res);

        policyPtr_->restrictResidual(res, coarseRes);

    }
    else
    {
        // No pre-sweeps done: x = 0 and residual = b
        policyPtr_->restrictResidual(b, coarseRes);
    }
}


template<class Type>
void Foam::coarseBlockAmgLevel<Type>::prolongateCorrection
(
    Field<Type>& x,
    const Field<Type>& coarseX
) const
{
    policyPtr_->prolongateCorrection(x, coarseX);
}


template<class Type>
void Foam::coarseBlockAmgLevel<Type>::smooth
(
    Field<Type>& x,
    const Field<Type>& b,
    const label nSweeps
) const
{
    smootherPtr_->smooth(x, b, nSweeps);
}


template<class Type>
void Foam::coarseBlockAmgLevel<Type>::solve
(
    Field<Type>& x,
    const Field<Type>& b,
    const scalar tolerance,
    const scalar relTol
) const
{

    if (matrixPtr_->symmetric())
    {
        BlockSolverPerformance<Type> coarseSolverPerf =
        //BlockGaussSeidelSolver<Type>
        BlockCGSolver<Type>
        (
            "topLevelCorr",
            matrixPtr_,
            dict_
        ).solve(x, b);

        if (lduMatrix::debug >= 2)
        {
            coarseSolverPerf.print();
        }
    }
    else
    {
        BlockSolverPerformance<Type> coarseSolverPerf =
        //BlockGaussSeidelSolver<Type>
        BlockBiCGStabSolver<Type>
        (
            "topLevelCorr",
            matrixPtr_,
            dict_
        ).solve(x, b);

        if (lduMatrix::debug >= 2)
        {
            coarseSolverPerf.print();
        }
    }

}


template<class Type>
void Foam::coarseBlockAmgLevel<Type>::scaleX
(
    Field<Type>& x,
    const Field<Type>& b,
    Field<Type>& xBuffer
) const
{

    // KRJ: 2013-02-05: Creating a new field not re-using
    Field<Type> Ax(x.size());

    matrixPtr_->Amul
    (
        reinterpret_cast<Field<Type>&>(Ax),
        x
    );

    scalar scalingFactorNum = sumProd(x,b);
    scalar scalingFactorDenom = sumProd(x,Ax);

    vector2D scalingVector(scalingFactorNum, scalingFactorDenom);
    reduce(scalingVector, sumOp<vector2D>());

    // Scale x
    if
    (
        mag(scalingVector[0]) > GREAT
     || mag(scalingVector[1]) > GREAT
     || scalingVector[0]*scalingVector[1] <= 0
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
Foam::autoPtr<Foam::BlockAmgLevel<Type> > Foam::coarseBlockAmgLevel<Type>::makeNextLevel() const
{
    if (policyPtr_->coarsen())
    {
        return autoPtr<Foam::BlockAmgLevel<Type> >
        (
            new coarseBlockAmgLevel
            (
                policyPtr_->restrictMatrix(),
                dict(),
                policyPtr_->type(),
                policyPtr_->groupSize(),
                policyPtr_->minCoarseEqns(),
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
