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
    BlockILUSolver

Description
    ILU solver

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    BlockILUSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BlockILUSolver_H
#define BlockILUSolver_H

#include "blockLduSolvers.H"
#include "BlockIterativeSolver.H"
#include "BlockSolverPerformance.H"
#include "blockCholeskyPrecons.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class BlockILUSolver Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockILUSolver
:
    public BlockIterativeSolver<Type>
{
    // Private data

        //- Gauss-Seidel preconditioner
        BlockCholeskyPrecon<Type> ilu_;

        //- Number of sweeps before evaluating residual
        label nSweeps_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockILUSolver(const BlockILUSolver<Type>&);

        //- Disallow default bitwise assignment
        void operator=(const BlockILUSolver<Type>&);


public:

    //- Runtime type information
    TypeName("ILU");


    // Constructors

        //- Construct from matrix components and solver data stream
        BlockILUSolver
        (
            const word& fieldName,
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        );


        //- Destructor
        virtual ~BlockILUSolver()
        {}


    // Member Functions

        //- Solve the matrix with this solver
        virtual BlockSolverPerformance<Type> solve
        (
            Field<Type>& x,
            const Field<Type>& b
        );


        //- Re-initialise solver after matrix coefficient update
        virtual void initMatrix()
        {
            ilu_.initMatrix();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BlockILUSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
