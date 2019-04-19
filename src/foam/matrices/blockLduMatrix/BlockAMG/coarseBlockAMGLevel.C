/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    coarseBlockAMGLevel

Description
    Coarse AMG level stores matrix, x and b locally, for BlockLduMatrix

Author
    Klas Jareteg, 2012-12-13

\*---------------------------------------------------------------------------*/

#include "coarseBlockAMGLevel.H"
#include "SubField.H"
#include "vector2D.H"
#include "coeffFields.H"
#include "BlockSolverPerformance.H"
#include "BlockCGSolver.H"
#include "BlockBiCGStabSolver.H"
#include "BlockGMRESSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Type>
Foam::coarseBlockAMGLevel<Type>::coarseBlockAMGLevel
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
    x_(matrixPtr_->diag().size(), pTraits<Type>::zero),
    b_(matrixPtr_->diag().size(), pTraits<Type>::zero),
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
//             ,"coarseSmoother"
        )
    ),
    Ax_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::coarseBlockAMGLevel<Type>::~coarseBlockAMGLevel()
{
    // Clear addressing interfaces
    if (addrPtr_.valid())
    {
        addrPtr_().clearInterfaces();
    }

    // Clear matrix interfaces
    if (matrixPtr_.valid())
    {
        matrixPtr_().clearInterfaces();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::BlockLduMatrix<Type>& Foam::coarseBlockAMGLevel<Type>::matrix()
{
    return matrixPtr_();
}


template<class Type>
Foam::Field<Type>& Foam::coarseBlockAMGLevel<Type>::x()
{
    return x_;
}


template<class Type>
Foam::Field<Type>& Foam::coarseBlockAMGLevel<Type>::b()
{
    return b_;
}


template<class Type>
void Foam::coarseBlockAMGLevel<Type>::residual
(
    const Field<Type>& x,
    const Field<Type>& b,
    Field<Type>& res
) const
{
    // Calculate residual
    matrixPtr_->Amul(res, x);

    // residual = b - Ax
    forAll (b, i)
    {
        res[i] = b[i] - res[i];
    }
}


template<class Type>
void Foam::coarseBlockAMGLevel<Type>::restrictResidual
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
void Foam::coarseBlockAMGLevel<Type>::prolongateCorrection
(
    Field<Type>& x,
    const Field<Type>& coarseX
) const
{
    coarseningPtr_->prolongateCorrection(x, coarseX);
}


template<class Type>
void Foam::coarseBlockAMGLevel<Type>::smooth
(
    Field<Type>& x,
    const Field<Type>& b,
    const label nSweeps
) const
{
    smootherPtr_->smooth(x, b, nSweeps);
}


template<class Type>
void Foam::coarseBlockAMGLevel<Type>::solve
(
    Field<Type>& x,
    const Field<Type>& b,
    const scalar tolerance,
    const scalar relTol
) const
{
    BlockSolverPerformance<Type> coarseSolverPerf
    (
        BlockBiCGStabSolver<Type>::typeName,
        "topLevelCorr"
    );

    // Create artificial dictionary for top-level solution
    dictionary topLevelDict;
    topLevelDict.add("nDirections", "5");
    topLevelDict.add("preconditioner", "ILUC0");
    topLevelDict.add("minIter", 0);
    topLevelDict.add("maxIter", 500);
    topLevelDict.add("tolerance", tolerance);
    topLevelDict.add("relTol", relTol);

    // Avoid issues with round-off on strict tolerance setup
    // HJ, 27/Jun/2013
    x = pTraits<Type>::zero;

    // Switch of debug in top-level direct solve
    label oldDebug = blockLduMatrix::debug();

    if (blockLduMatrix::debug >= 4)
    {
        blockLduMatrix::debug = 2;
    }
    else if (blockLduMatrix::debug == 3)
    {
        blockLduMatrix::debug = 1;
    }
    else
    {
        blockLduMatrix::debug = 0;
    }

    if (matrixPtr_->symmetric())
    {
        coarseSolverPerf = BlockCGSolver<Type>
        (
            "topLevelCorr",
            matrixPtr_,
            topLevelDict
        ).solve(x, b);
    }
    else
    {
        coarseSolverPerf =
        BlockBiCGStabSolver<Type>
        (
            "topLevelCorr",
            matrixPtr_,
            topLevelDict
        ).solve(x, b);
    }

    // Check for convergence

//    blockLduMatrix::debug = 1;
//    coarseSolverPerf.print();
//    Info<< "b: " << gSum(b) << " corr (min, max): (" << gMin(x) << ", " << gMax(x) << endl;

    const scalar magInitialRes = mag(coarseSolverPerf.initialResidual());
    const scalar magFinalRes = mag(coarseSolverPerf.finalResidual());

    if (magFinalRes > magInitialRes && magInitialRes > 1e-12)
    {
        if (blockLduMatrix::debug)
        {
            Info<< "Divergence in top AMG level" << endl;
            coarseSolverPerf.print();
        }

        x = pTraits<Type>::zero;
    }

    // Restore debug
    blockLduMatrix::debug = oldDebug;

    if (blockLduMatrix::debug >= 3)
    {
        coarseSolverPerf.print();
    }
}


template<class Type>
void Foam::coarseBlockAMGLevel<Type>::scaleX
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

    // Variant 1: scale complete x with a single scaling factor
    scalar scalingFactorNum = sumProd(x, b);
    scalar scalingFactorDenom = sumProd(x, Ax_);

    vector2D scalingVector(scalingFactorNum, scalingFactorDenom);
    reduce(scalingVector, sumOp<vector2D>());

    // Scale x
    if
    (
        mag(scalingVector[0]) > SMALL
     && mag(scalingVector[0]) < GREAT
     && mag(scalingVector[0]) > SMALL
     && mag(scalingVector[1]) < GREAT
     && scalingVector[0]*scalingVector[1] > 0
    )
    {
        // Regular scaling
        x *= Foam::max
        (
            0.01,
            Foam::min
            (
                scalingVector[0]/scalingVector[1],
                100
            )
        );
    }
}


template<class Type>
Foam::autoPtr<Foam::BlockAMGLevel<Type> >
Foam::coarseBlockAMGLevel<Type>::makeNextLevel() const
{
    if (coarseningPtr_->coarsen())
    {
        return coarseningPtr_->restrictMatrix();
    }
    else
    {
        // Final level: cannot coarsen
        return autoPtr<BlockAMGLevel<Type> >();
    }
}


template<class Type>
void Foam::coarseBlockAMGLevel<Type>::initLevel
(
    autoPtr<Foam::BlockAMGLevel<Type> >& coarseLevelPtr
)
{
    // Update matrix by repeating coarsening

    // Update smoother for new matrix
    smootherPtr_->initMatrix();

    // Update coarse matrix if it exists
    if (coarseLevelPtr.valid())
    {
        coarseningPtr_->updateMatrix(coarseLevelPtr->matrix());
    }
}


// ************************************************************************* //
