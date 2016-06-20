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

Class
    fineBlockAMGLevel

Description
    Finest AMG level container, using matrix and field references

Author
    Klas Jareteg, 2013-04-15

\*---------------------------------------------------------------------------*/

#include "fineBlockAMGLevel.H"
#include "coarseBlockAMGLevel.H"
#include "BlockSolverPerformance.H"
#include "BlockCoeffNorm.H"
#include "BlockCoeffTwoNorm.H"
#include "BlockBiCGStabSolver.H"
#include "BlockCGSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Type>
Foam::fineBlockAMGLevel<Type>::fineBlockAMGLevel
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict,
    const word& coarseningType,
    const label groupSize,
    const label minCoarseEqns
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
//             ,"fineSmoother"
        )
    ),
    Ax_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::Field<Type>& Foam::fineBlockAMGLevel<Type>::x()
{
    FatalErrorIn("Field<Type>& Foam::fineBlockAMGLevel<Type>::x()")
        << "x is not available."
        << abort(FatalError);

    // Dummy return
    return const_cast<Field<Type>&>(Field<Type>::null());
}


template<class Type>
Foam::Field<Type>& Foam::fineBlockAMGLevel<Type>::b()
{
    FatalErrorIn("Field<Type>& Foam::fineBlockAMGLevel<Type>::b()")
        << "b is not available."
        << abort(FatalError);

    // Dummy return
    return const_cast<Field<Type>&>(Field<Type>::null());
}


template<class Type>
void Foam::fineBlockAMGLevel<Type>::residual
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
void Foam::fineBlockAMGLevel<Type>::restrictResidual
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
void Foam::fineBlockAMGLevel<Type>::prolongateCorrection
(
    Field<Type>& x,
    const Field<Type>& coarseX
) const
{
    coarseningPtr_->prolongateCorrection(x, coarseX);
}


template<class Type>
void Foam::fineBlockAMGLevel<Type>::smooth
(
    Field<Type>& x,
    const Field<Type>& b,
    const label nSweeps
) const
{
    smootherPtr_->smooth(x, b, nSweeps);
}


template<class Type>
void Foam::fineBlockAMGLevel<Type>::solve
(
    Field<Type>& x,
    const Field<Type>& b,
    const scalar tolerance,
    const scalar relTol
) const
{
    Info<< "Fine level solver" << endl;

    // Create artificial dictionary for finest solution
    dictionary finestDict;
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
    }
}


template<class Type>
void Foam::fineBlockAMGLevel<Type>::scaleX
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

    matrix_.Amul(Ax_, x);

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
Foam::fineBlockAMGLevel<Type>::makeNextLevel() const
{
    if (coarseningPtr_->coarsen())
    {
        return coarseningPtr_->restrictMatrix();
    }
    else
    {
        // Final level: cannot coarsen
        return autoPtr<Foam::BlockAMGLevel<Type> >();
    }
}


// ************************************************************************* //
