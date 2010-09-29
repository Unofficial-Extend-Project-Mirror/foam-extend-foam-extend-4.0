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
