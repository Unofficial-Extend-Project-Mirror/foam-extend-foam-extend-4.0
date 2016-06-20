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

\*---------------------------------------------------------------------------*/

#include "FDICPreconditioner.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FDICPreconditioner, 0);

    lduMatrix::preconditioner::
        addsymMatrixConstructorToTable<FDICPreconditioner>
        addFDICPreconditionerSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FDICPreconditioner::FDICPreconditioner
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary&
)
:
    lduPreconditioner
    (
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces
    ),
    rD_(matrix.diag()),
    rDuUpper_(matrix.upper().size()),
    rDlUpper_(matrix.upper().size())
{
    scalar* __restrict__ rDPtr = rD_.begin();
    scalar* __restrict__ rDuUpperPtr = rDuUpper_.begin();
    scalar* __restrict__ rDlUpperPtr = rDlUpper_.begin();

    const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        matrix_.lduAddr().lowerAddr().begin();
    const scalar* const __restrict__ upperPtr = matrix_.upper().begin();

    register label nCells = rD_.size();
    register label nFaces = matrix_.upper().size();

    for (register label face=0; face<nFaces; face++)
    {
        rDPtr[uPtr[face]] -= sqr(upperPtr[face])/rDPtr[lPtr[face]];
    }

    // Generate reciprocal FDIC
    for (register label cell=0; cell<nCells; cell++)
    {
        rDPtr[cell] = 1.0/rDPtr[cell];
    }

    for (register label face=0; face<nFaces; face++)
    {
        rDuUpperPtr[face] = rDPtr[uPtr[face]]*upperPtr[face];
        rDlUpperPtr[face] = rDPtr[lPtr[face]]*upperPtr[face];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FDICPreconditioner::precondition
(
    scalarField& wA,
    const scalarField& rA,
    const direction
) const
{
    scalar* __restrict__ wAPtr = wA.begin();
    const scalar* __restrict__ rAPtr = rA.begin();
    const scalar* __restrict__ rDPtr = rD_.begin();

    const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        matrix_.lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ rDuUpperPtr = rDuUpper_.begin();
    const scalar* const __restrict__ rDlUpperPtr = rDlUpper_.begin();

    register label nCells = wA.size();
    register label nFaces = matrix_.upper().size();
    register label nFacesM1 = nFaces - 1;

    for (register label cell=0; cell<nCells; cell++)
    {
        wAPtr[cell] = rDPtr[cell]*rAPtr[cell];
    }

    for (register label face=0; face<nFaces; face++)
    {
        wAPtr[uPtr[face]] -= rDuUpperPtr[face]*wAPtr[lPtr[face]];
    }

    for (register label face=nFacesM1; face>=0; face--)
    {
        wAPtr[lPtr[face]] -= rDlUpperPtr[face]*wAPtr[uPtr[face]];
    }
}


// ************************************************************************* //
