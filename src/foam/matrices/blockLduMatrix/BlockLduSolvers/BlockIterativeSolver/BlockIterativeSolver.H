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
    BlockIterativeSolver

Description
    Virtual base class for a block iterative solvers

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    BlockIterativeSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BlockIterativeSolver_H
#define BlockIterativeSolver_H

#include "BlockLduSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class BlockIterativeSolver Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockIterativeSolver
:
    public BlockLduSolver<Type>
{
    // Private data

        //- Solver tolerance
        scalar tolerance_;

        //- Relative tolerance
        scalar relTolerance_;

        //- Minimum number of iterations
        //  (forced irrespective of convergence)
        label minIter_;

        //- Maximum number of iterations
        label maxIter_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockIterativeSolver(const BlockIterativeSolver&);

        //- Disallow default bitwise assignment
        void operator=(const BlockIterativeSolver&);


protected:

    // Protected Member Functions

        //- Return normalisation factor
        scalar normFactor
        (
            Field<Type>& x,
            const Field<Type>& b
        ) const;

        //- Is the stop criterion reached
        bool stop(BlockSolverPerformance<Type>& solverPerf) const;


public:

    // Constructors

        //- Construct from matrix components and solver data stream
        BlockIterativeSolver
        (
            const word& fieldName,
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        );


    // Destructor

        virtual ~BlockIterativeSolver()
        {}


    // Member Functions

        // Access

            //- Return tolerance
            scalar tolerance() const
            {
                return tolerance_;
            }

            //- Return relative tolerance
            scalar relTolerance() const
            {
                return relTolerance_;
            }

            //- Return minimum number of iterations
            label minIter() const
            {
                return minIter_;
            }

            //- Return maximum number of iterations
            label maxIter() const
            {
                return maxIter_;
            }


        // Solve

            //- Solve the matrix with this solver
            virtual BlockSolverPerformance<Type> solve
            (
                Field<Type>& x,
                const Field<Type>& b
            ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "BlockIterativeSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
