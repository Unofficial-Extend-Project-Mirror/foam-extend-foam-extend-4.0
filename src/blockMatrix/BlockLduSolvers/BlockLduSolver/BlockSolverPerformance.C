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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "BlockSolverPerformance.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::BlockSolverPerformance<Type>::checkConvergence
(
    const scalar Tolerance,
    const scalar RelTolerance
)
{
    if (BlockLduMatrix<Type>::debug >= 2)
    {
        Info<< solverName_
            << ":  Iteration " << nIterations_
            << " residual = " << finalResidual_
            << endl;
    }

    if
    (
        mag(finalResidual_) < Tolerance
     || (
            RelTolerance > SMALL
         && mag(finalResidual_) <= RelTolerance*mag(initialResidual_)
        )
    )
    {
        converged_ = true;
    }
    else
    {
        converged_ = false;
    }

    return converged_;
}


template<class Type>
bool Foam::BlockSolverPerformance<Type>::checkSingularity
(
    const scalar& residual
)
{
    if (mag(residual) > VSMALL)
    {
        singular_ = false;
    }
    else
    {
        singular_ = true;
    }

    return singular_;
}


template<class Type>
void Foam::BlockSolverPerformance<Type>::print() const
{
    if (BlockLduMatrix<Type>::debug)
    {
        Info<< solverName_ << ":  Solving for " << fieldName_;

        if (singular())
        {
            Info<< ":  solution singularity" << endl;
        }
        else
        {
            Info<< ", Initial residual = " << initialResidual_
                << ", Final residual = " << finalResidual_
                << ", No Iterations " << nIterations_
                << endl;
        }
    }
}


// ************************************************************************* //
