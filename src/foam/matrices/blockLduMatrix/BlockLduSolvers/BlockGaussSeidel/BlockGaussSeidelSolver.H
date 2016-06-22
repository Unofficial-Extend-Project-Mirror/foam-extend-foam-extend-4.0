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
    BlockGaussSeidelSolver

Description
    Gauss-Seidel solver

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    BlockGaussSeidelSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BlockGaussSeidelSolver_H
#define BlockGaussSeidelSolver_H

#include "blockLduSolvers.H"
#include "BlockIterativeSolver.H"
#include "BlockSolverPerformance.H"
#include "blockGaussSeidelPrecons.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class BlockGaussSeidelSolver Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockGaussSeidelSolver
:
    public BlockIterativeSolver<Type>
{
    // Private data

        //- Gauss-Seidel preconditioner
        BlockGaussSeidelPrecon<Type> gs_;

        //- Number of sweeps before evaluating residual
        label nSweeps_;



    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockGaussSeidelSolver(const BlockGaussSeidelSolver<Type>&);

        //- Disallow default bitwise assignment
        void operator=(const BlockGaussSeidelSolver<Type>&);


public:

    //- Runtime type information
    TypeName("GaussSeidel");


    // Constructors

        //- Construct from matrix components and solver data stream
        BlockGaussSeidelSolver
        (
            const word& fieldName,
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        );


        // Destructor

            virtual ~BlockGaussSeidelSolver()
            {}


    // Member Functions

        //- Solve the matrix with this solver
        virtual BlockSolverPerformance<Type> solve
        (
            Field<Type>& x,
            const Field<Type>& b
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockGaussSeidelSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
