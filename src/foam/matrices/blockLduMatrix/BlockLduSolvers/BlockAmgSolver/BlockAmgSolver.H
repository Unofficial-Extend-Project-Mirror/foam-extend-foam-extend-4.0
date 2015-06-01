/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    BlockAmgSolver

Description
    Algebraic Multigrid solver with run-time selection of coarsening and cycle
    adjusted for BlockLduMatrix.

Author
    Klas Jareteg, 2012-12-13

SourceFiles
    BlockAmgSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BlockAmgSolver_H
#define BlockAmgSolver_H

#include "BlockAmgPrecon.H"
#include "blockLduSolvers.H"
#include "blockLduPrecons.H"
#include "BlockIterativeSolver.H"
#include "BlockSolverPerformance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class BlockAmgSolver Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class BlockAmgSolver
:
    public BlockIterativeSolver<Type>
{
private:

    // Private data

        //- AMG preconditioner
        BlockAmgPrecon<Type> amg_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockAmgSolver(const BlockAmgSolver<Type>&);

        //- Disallow default bitwise assignment
        void operator=(const BlockAmgSolver<Type>&);


public:

    //- Runtime type information
    TypeName("AMG");


    // Constructors

        //- Construct from matrix components and solver data stream
        BlockAmgSolver
        (
            const word& fieldName,
            const BlockLduMatrix<Type>& matrix,
            const dictionary& dict
        );


    // Destructor

        virtual ~BlockAmgSolver()
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
#   include "BlockAmgSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
