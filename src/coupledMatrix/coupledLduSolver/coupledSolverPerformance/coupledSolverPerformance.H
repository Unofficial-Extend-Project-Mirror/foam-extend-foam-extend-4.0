/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    coupledSolverPerformance

Description
    Class returned by a coupled solver.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    coupledSolverPerformance.C

\*---------------------------------------------------------------------------*/

#ifndef coupledSolverPerformance_H
#define coupledSolverPerformance_H

#include "coupledLduMatrix.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class coupledSolverPerformance Declaration
\*---------------------------------------------------------------------------*/

class coupledSolverPerformance
{
    // Private data

        //- Solver name
        word solverName_;

        //- Field name
        word fieldName_;

        //- Initial residual
        scalar initialResidual_;

        //- Final residual
        scalar finalResidual_;

        //- Number of iterations
        label nIterations_;

        //- Has the solver converged?
        bool converged_;

        //- Is the matrix singular?
        bool singular_;


public:

    // Constructors

        //- Construct with solver and field name
        coupledSolverPerformance
        (
            const word& solverName,
            const word& fieldName,
            const scalar iRes = 0,
            const scalar fRes = 0,
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
        const scalar& initialResidual() const
        {
            return initialResidual_;
        }

        //- Return initial residual
        scalar& initialResidual()
        {
            return initialResidual_;
        }

        //- Return final residual
        const scalar& finalResidual() const
        {
            return finalResidual_;
        }

        //- Return final residual
        scalar& finalResidual()
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

#endif

// ************************************************************************* //
