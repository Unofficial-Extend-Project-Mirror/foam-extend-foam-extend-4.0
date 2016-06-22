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
    BlockAMGCycle

Description
    Algebraic multigrid cycle class for BlockLduMatrix

Author
    Klas Jareteg, 2012-12-12

\*---------------------------------------------------------------------------*/

#include "BlockAMGCycle.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from AMG level
template<class Type>
Foam::BlockAMGCycle<Type>::BlockAMGCycle
(
    autoPtr<BlockAMGLevel<Type> > levelPtr
)
:
    levelPtr_(levelPtr),
    coarseLevelPtr_(NULL),
    nLevels_(0)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockAMGCycle<Type>::~BlockAMGCycle()
{
    deleteDemandDrivenData(coarseLevelPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::BlockAMGCycle<Type>::makeCoarseLevels(const label nMaxLevels)
{
    // Make coarse levels
    if (nLevels_ == 0)
    {
        bool addCoarse = true;
        BlockAMGCycle<Type>* curCyclePtr = this;

        // Do forever
        for (;;)
        {
            nLevels_++;

            autoPtr<BlockAMGLevel<Type> > coarsePtr =
                curCyclePtr->levelPtr_->makeNextLevel();

            // Check if a coarse level is valid and allowed
            if (nLevels_ >= nMaxLevels || !coarsePtr.valid())
            {
                addCoarse = false;
            }

            reduce(addCoarse, andOp<bool>());

            if (addCoarse)
            {
                curCyclePtr->coarseLevelPtr_ =
                    new BlockAMGCycle<Type>(coarsePtr);

                // Point to the next coarse level
                curCyclePtr = curCyclePtr->coarseLevelPtr_;
            }
            else
            {
                break;
            }
        }

        if (BlockLduMatrix<Type>::debug >= 2)
        {
            Info<< "Created " << nLevels_ << " AMG levels" << endl;
        }
    }
}


template<class Type>
void Foam::BlockAMGCycle<Type>::fixedCycle
(
    Field<Type>& x,
    const Field<Type>& b,
    Field<Type>& xBuffer,
    const blockAMGCycleName::cycleType cycle,
    const label nPreSweeps,
    const label nPostSweeps,
    const bool scale
) const
{
    if (coarseLevelPtr_)
    {
        // Pre-smoothing
        levelPtr_->smooth(x, b, nPreSweeps);

        // Get reference to coarse level
        Field<Type>& xCoarse = coarseLevelPtr_->levelPtr_->x();
        Field<Type>& bCoarse = coarseLevelPtr_->levelPtr_->b();

        // Zero out coarse x
        xCoarse = pTraits<Type>::zero;

        // Restrict residual: optimisation on number of pre-sweeps
        levelPtr_->restrictResidual
        (
            x,
            b,
            xBuffer,
            bCoarse,
            nPreSweeps > 0 || cycle != V_CYCLE
        );

        coarseLevelPtr_->fixedCycle
        (
            xCoarse,
            bCoarse,
            xBuffer,
            cycle,
            nPreSweeps,
            nPostSweeps,
            scale
        );

        if (cycle == F_CYCLE)
        {
            coarseLevelPtr_->fixedCycle
            (
                xCoarse,
                bCoarse,
                xBuffer,
                V_CYCLE,
                nPreSweeps,
                nPostSweeps,
                scale
            );
        }
        else if (cycle == W_CYCLE)
        {
            coarseLevelPtr_->fixedCycle
            (
                xCoarse,
                bCoarse,
                xBuffer,
                W_CYCLE,
                nPreSweeps,
                nPostSweeps,
                scale
            );
        }

        if (scale)
        {
            // Calculate scaling factor using a buffer
            coarseLevelPtr_->levelPtr_->scaleX(xCoarse, bCoarse, xBuffer);
        }

        levelPtr_->prolongateCorrection(x, xCoarse);

        // Post-smoothing
        levelPtr_->smooth(x, b, nPostSweeps);
    }
    else
    {
        // Call direct solver
        levelPtr_->solve(x, b, 1e-9, 0);
    }
}


// ************************************************************************* //
