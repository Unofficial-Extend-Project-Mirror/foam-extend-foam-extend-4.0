/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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

Class
    BlockAmgCycle

Description
    Algebraic multigrid cycle class for BlockLduMatrix

Author
    Klas Jareteg, 2012-12-12

\*---------------------------------------------------------------------------*/

#include "BlockAmgCycle.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from AMG level
template<class Type>
Foam::BlockAmgCycle<Type>::BlockAmgCycle
(
    autoPtr<BlockAmgLevel<Type> > levelPtr
)
:
    levelPtr_(levelPtr),
    coarseLevelPtr_(NULL),
    nLevels_(0)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockAmgCycle<Type>::~BlockAmgCycle()
{
    deleteDemandDrivenData(coarseLevelPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::BlockAmgCycle<Type>::makeCoarseLevels(const label nMaxLevels)
{
    // Make coarse levels
    if (nLevels_ == 0)
    {
        bool addCoarse = true;
        BlockAmgCycle<Type>* curCyclePtr = this;

        // Do forever
        for (;;)
        {
            nLevels_++;

            autoPtr<BlockAmgLevel<Type> > coarsePtr =
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
                    new BlockAmgCycle<Type>(coarsePtr);

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
void Foam::BlockAmgCycle<Type>::fixedCycle
(
    Field<Type>& x,
    const Field<Type>& b,
    Field<Type>& xBuffer,
    const blockAmgCycleName::cycleType cycle,
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
        forAll (xCoarse,i)
        {
            xCoarse[i] *= 0.0;
        }

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
