/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Class
    coarseBlockAmgLevel

Description
    Coarse AMG level stores matrix, x and b locally, for BlockLduMatrix

Author
    Klas Jareteg, 2012-12-13

\*---------------------------------------------------------------------------*/

#include "coarseBlockAmgLevel.H"
#include "SubField.H"
#include "vector2D.H"
#include "coeffFields.H"
#include "BlockSolverPerformance.H"
// #include "BlockBiCGStabSolver.H"
// #include "BlockCGSolver.H"
#include "BlockGMRESSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Type>
Foam::coarseBlockAmgLevel<Type>::coarseBlockAmgLevel
(
    autoPtr<lduPrimitiveMesh> addrPtr,
    autoPtr<BlockLduMatrix<Type> > matrixPtr,
    const dictionary& dict,
    const word& coarseningType,
    const label groupSize,
    const label minCoarseEqns
)
:
    addrPtr_(addrPtr),
    matrixPtr_(matrixPtr),
    x_(matrixPtr_->diag().size(),pTraits<Type>::zero),
    b_(matrixPtr_->diag().size(),pTraits<Type>::zero),
    dict_(dict),
    coarseningPtr_
    (
        BlockMatrixCoarsening<Type>::New
        (
            coarseningType,
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
    ),
    Ax_()
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

        coarseningPtr_->restrictResidual(res, coarseRes);
    }
    else
    {
        // No pre-sweeps done: x = 0 and residual = b
        coarseningPtr_->restrictResidual(b, coarseRes);
    }
}


template<class Type>
void Foam::coarseBlockAmgLevel<Type>::prolongateCorrection
(
    Field<Type>& x,
    const Field<Type>& coarseX
) const
{
    coarseningPtr_->prolongateCorrection(x, coarseX);
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
    BlockSolverPerformance<Type> coarseSolverPerf
    (
        BlockGMRESSolver<Type>::typeName,
        "topLevelCorr"
    );

    label maxIter = Foam::min(2*coarseningPtr_->minCoarseEqns(), 1000);

    // Create artificial dictionary for top-level solution
    dictionary topLevelDict;
    topLevelDict.add("nDirections", "5");
    topLevelDict.add("minIter", 1);
    topLevelDict.add("maxIter", maxIter);
    topLevelDict.add("tolerance", tolerance);
    topLevelDict.add("relTol", relTol);

    // Avoid issues with round-off on strict tolerance setup
    // HJ, 27/Jun/2013
    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

//     x = inverse(diag) & b
    CoeffField<Type> invDiag = inv(matrixPtr_->diag());
    multiply(x, invDiag, b);

    // Do not solve if the number of equations is smaller than 5
    if (coarseningPtr_->minCoarseEqns() < 5)
    {
        return;
    }

    // Switch of debug in top-level direct solve
    label oldDebug = BlockLduMatrix<Type>::debug();

    if (BlockLduMatrix<Type>::debug >= 4)
    {
        BlockLduMatrix<Type>::debug = 1;
    }
    else
    {
        BlockLduMatrix<Type>::debug = 0;
    }

    if (matrixPtr_->symmetric())
    {
        topLevelDict.add("preconditioner", "Cholesky");

        coarseSolverPerf =
//         BlockCGSolver<Type>
        BlockGMRESSolver<Type>
        (
            "topLevelCorr",
            matrixPtr_,
            topLevelDict
        ).solve(x, b);
    }
    else
    {
        topLevelDict.add("preconditioner", "Cholesky");

        coarseSolverPerf =
//         BlockBiCGStabSolver<Type>
        BlockGMRESSolver<Type>
        (
            "topLevelCorr",
            matrixPtr_,
            topLevelDict
        ).solve(x, b);
    }

    // Restore debug
    BlockLduMatrix<Type>::debug = oldDebug;

    // Escape cases of top-level solver divergence
    if
    (
        coarseSolverPerf.nIterations() == maxIter
     && (
            coarseSolverPerf.finalResidual()
         >= coarseSolverPerf.initialResidual()
        )
    )
    {
        // Top-level solution failed.  Attempt rescue
        // HJ, 27/Jul/2013
        multiply(x, invDiag, b);

        // Print top level correction failure as info for user
        coarseSolverPerf.print();
    }

    if (BlockLduMatrix<Type>::debug >= 2)
    {
        coarseSolverPerf.print();
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
    // KRJ: 2013-02-05: Removed subfield, creating a new field
    // Initialise and size buffer to avoid multiple allocation.
    // Buffer is created as private data of AMG level
    // HJ, 26/Feb/2015
    if (Ax_.empty())
    {
        Ax_.setSize(x.size());
    }

    matrixPtr_->Amul(Ax_, x);

#if 0

    // Variant 1 (old): scale complete x with a single scaling factor
    scalar scalingFactorNum = sumProd(x, b);
    scalar scalingFactorDenom = sumProd(x, Ax_);

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
        x *= scalingVector[0]/stabilise(scalingVector[1], VSMALL);
    }

#else

    // Variant 2: scale each x individually
    // HJ, 25/Feb/2015
    Type scalingFactorNum = sumCmptProd(x, b);
    Type scalingFactorDenom = sumCmptProd(x, Ax_);

    Vector2D<Type> scalingVector(scalingFactorNum, scalingFactorDenom);
    reduce(scalingVector, sumOp<Vector2D<Type> >());

    Type scalingFactor = pTraits<Type>::one;

    // Scale x
    for (direction dir = 0; dir < pTraits<Type>::nComponents; dir++)
    {
        scalar num = component(scalingVector[0], dir);
        scalar denom = component(scalingVector[1], dir);

        if
        (
            mag(num) > GREAT || mag(denom) > GREAT
         || num*denom <= 0 || mag(num) < mag(denom)
        )
        {
            // Factor = 1.0, no scaling
        }
        else if (mag(num) > 2*mag(denom))
        {
            setComponent(scalingFactor, dir) = 2;
        }
        else
        {
            // Regular scaling
            setComponent(scalingFactor, dir) = num/stabilise(denom, VSMALL);
        }
    }

    // Scale X
    cmptMultiply(x, x, scalingFactor);

#endif
}


template<class Type>
Foam::autoPtr<Foam::BlockAmgLevel<Type> >
Foam::coarseBlockAmgLevel<Type>::makeNextLevel() const
{
    if (coarseningPtr_->coarsen())
    {
        return coarseningPtr_->restrictMatrix();
    }
    else
    {
        // Final level: cannot coarsen
        return autoPtr<BlockAmgLevel<Type> >();
    }
}


// ************************************************************************* //
