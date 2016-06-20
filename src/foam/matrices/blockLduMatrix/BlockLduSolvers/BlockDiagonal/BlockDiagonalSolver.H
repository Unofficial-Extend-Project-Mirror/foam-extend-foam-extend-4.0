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
    BlockDiagonalSolver

Description
    Solver for diagonal matrices

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    BlockDiagonalSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BlockDiagonalSolver_H
#define BlockDiagonalSolver_H

#include "Field.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class BlockLduMatrix;

template<class Type>
class BlockLduSolver;

template<class Type>
class BlockSolverPerformance;

/*---------------------------------------------------------------------------*\
                      Class BlockDiagonalSolver Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockDiagonalSolver
:
    public BlockLduSolver<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockDiagonalSolver(const BlockDiagonalSolver<Type>&);

        //- Disallow default bitwise assignment
        void operator=(const BlockDiagonalSolver<Type>&);


public:

    //- Runtime type information
    TypeName("Diagonal");


    // Constructors

        //- Construct from matrix
        BlockDiagonalSolver
        (
            const word& fieldName,
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        );


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
#   include "BlockDiagonalSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
