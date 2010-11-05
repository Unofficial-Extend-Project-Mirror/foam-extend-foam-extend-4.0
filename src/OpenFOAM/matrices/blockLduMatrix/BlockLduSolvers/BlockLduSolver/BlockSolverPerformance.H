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

Class
    BlockSolverPerformance

Description
    Class returned by a solver.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    BlockSolverPerformance.C

\*---------------------------------------------------------------------------*/

#ifndef BlockSolverPerformance_H
#define BlockSolverPerformance_H

#include "error.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //

template<class Type>
class BlockSolverPerformance;

template<class Type>
Ostream& operator<<(Ostream&, const BlockSolverPerformance<Type>&);


/*---------------------------------------------------------------------------*\
                   Class BlockSolverPerformance Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockSolverPerformance
{
    // Private data

        //- Solver name
        word solverName_;

        //- Field name
        word fieldName_;

        //- Initial residual
        Type initialResidual_;

        //- Final residual
        Type finalResidual_;

        //- Number of iterations
        label nIterations_;

        //- Has the solver converged?
        bool converged_;

        //- Is the matrix singular?
        bool singular_;


public:

    // Constructors

        //- Construct with solver and field name
        BlockSolverPerformance
        (
            const word& solverName,
            const word& fieldName,
            const Type& iRes = pTraits<Type>::zero,
            const Type& fRes = pTraits<Type>::zero,
            const label nIter = 0,
            const bool converged = false,
            const bool singular = false
        )
        :
            solverName_(solverName),
            fieldName_(fieldName),
            initialResidual_(iRes),
            finalResidual_(fRes),
            nIterations_(nIter),
            converged_(converged),
            singular_(singular)
        {}


    // Member functions

        //- Return solver name
        const word& solverName() const
        {
            return solverName_;
        }

        //- Return initial residual
        const Type& initialResidual() const
        {
            return initialResidual_;
        }

        //- Return initial residual
        Type& initialResidual()
        {
            return initialResidual_;
        }

        //- Return final residual
        const Type& finalResidual() const
        {
            return finalResidual_;
        }

        //- Return final residual
        Type& finalResidual()
        {
            return finalResidual_;
        }

        //- Return number of iterations
        label nIterations() const
        {
            return nIterations_;
        }

        //- Return number of iterations
        label& nIterations()
        {
            return nIterations_;
        }

        //- Has the solver converged?
        bool converged() const
        {
            return converged_;
        }

        //- Has the solver converged?
        bool& converged()
        {
            return converged_;
        }

        //- Is the matrix singular?
        bool singular() const
        {
            return singular_;
        }

        //- Is the matrix singular?
        bool& singular()
        {
            return singular_;
        }

        //- Convergence test
        bool checkConvergence
        (
            const scalar tolerance,
            const scalar relTolerance
        );

        //- Singularity test
        bool checkSingularity(const scalar& residual);

        //- Print summary of solver performance
        void print() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockSolverPerformance.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
