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
    amgCycle

Description
    Algebraic multigrid cycle class

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "amgCycle.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* Foam::NamedEnum<Foam::amgCycle::cycleType, 3>::names[] =
{
    "V-cycle",
    "W-cycle",
    "F-cycle"
};


const Foam::NamedEnum<Foam::amgCycle::cycleType, 3>
Foam::amgCycle::cycleNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from AMG level
Foam::amgCycle::amgCycle(autoPtr<amgLevel> levelPtr)
:
    levelPtr_(levelPtr),
    coarseLevelPtr_(NULL),
    nLevels_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::amgCycle::~amgCycle()
{
    deleteDemandDrivenData(coarseLevelPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::amgCycle::makeCoarseLevels(const label nMaxLevels)
{
    // Make coarse levels
    if (nLevels_ == 0)
    {
        bool addCoarse = true;
        amgCycle* curCyclePtr = this;

        // Do forever
        for (;;)
        {
            nLevels_++;

            autoPtr<amgLevel> coarsePtr =
                curCyclePtr->levelPtr_->makeNextLevel();

            // Check if a coarse level is valid and allowed
            if (nLevels_ >= nMaxLevels || !coarsePtr.valid())
            {
                addCoarse = false;
            }

            reduce(addCoarse, andOp<bool>());

            if (addCoarse)
            {
                curCyclePtr->coarseLevelPtr_ = new amgCycle(coarsePtr);

                // Point to the next coarse level
                curCyclePtr = curCyclePtr->coarseLevelPtr_;
            }
            else
            {
                break;
            }
        }

        if (lduMatrix::debug >= 2)
        {
            Info<< "Created " << nLevels_ << " AMG levels" << endl;
        }
    }
}


void Foam::amgCycle::fixedCycle
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt,
    scalarField& xBuffer,
    const cycleType cycle,
    const label nPreSweeps,
    const label nPostSweeps,
    const bool scale
) const
{
    if (coarseLevelPtr_)
    {
        // Pre-smoothing
        levelPtr_->smooth(x, b, cmpt, nPreSweeps);

        // Get reference to coarse level
        scalarField& xCoarse = coarseLevelPtr_->levelPtr_->x();
        scalarField& bCoarse = coarseLevelPtr_->levelPtr_->b();

        // Zero out coarse x
        xCoarse = 0;

        // Restrict residual: optimisation on number of pre-sweeps
        levelPtr_->restrictResidual
        (
            x,
            b,
            cmpt,
            xBuffer,
            bCoarse,
            nPreSweeps > 0 || cycle != V_CYCLE
        );

        coarseLevelPtr_->fixedCycle
        (
            xCoarse,
            bCoarse,
            cmpt,
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
                cmpt,
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
                cmpt,
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

            coarseLevelPtr_->levelPtr_->scaleX(xCoarse, bCoarse, cmpt, xBuffer);
        }

        levelPtr_->prolongateCorrection(x, xCoarse);

        // Post-smoothing
        levelPtr_->smooth(x, b, cmpt, nPostSweeps);
    }
    else
    {
        // Call direct solver
        levelPtr_->solve(x, b, cmpt, 1e-9, 0);
    }
}


// ************************************************************************* //
