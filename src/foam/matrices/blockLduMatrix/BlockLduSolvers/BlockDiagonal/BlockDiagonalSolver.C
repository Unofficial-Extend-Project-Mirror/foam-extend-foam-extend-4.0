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

Description
    Solver for diagonal matrices.

\*---------------------------------------------------------------------------*/

#include "BlockDiagonalSolver.H"
#include "BlockSolverPerformance.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
template<class Type>
Foam::BlockDiagonalSolver<Type>::BlockDiagonalSolver
(
    const word& fieldName,
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
:
    BlockLduSolver<Type>(fieldName, matrix, dict)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockSolverPerformance<Type>
Foam::BlockDiagonalSolver<Type>::solve
(
    Field<Type>& x,
    const Field<Type>& b
)
{
    CoeffField<Type> dD = inv(this->matrix_.diag());

    multiply(x, b, dD);

    return BlockSolverPerformance<Type>
    (
        this->typeName,
        this->fieldName(),
        pTraits<Type>::zero,
        pTraits<Type>::zero,
        0,
        true,
        false
    );
}


// ************************************************************************* //
