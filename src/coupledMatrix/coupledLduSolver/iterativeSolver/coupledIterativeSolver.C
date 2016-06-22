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
    Virtual base class for coupled iterative solvers

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "coupledIterativeSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledIterativeSolver::coupledIterativeSolver
(
    const word& fieldName,
    const coupledLduMatrix& matrix,
    const PtrList<FieldField<Field, scalar> >& bouCoeffs,
    const PtrList<FieldField<Field, scalar> >& intCoeffs,
    const lduInterfaceFieldPtrsListList& interfaces,
    const dictionary& solverData
)
:
    coupledLduSolver
    (
        fieldName,
        matrix,
        bouCoeffs,
        intCoeffs,
        interfaces
    ),
    dict_(solverData),
    tolerance_(1e-6),
    relTolerance_(0),
    minIter_(0),
    maxIter_(1000)
{
    readControls();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::coupledIterativeSolver::dict() const
{
    return dict_;
}


void Foam::coupledIterativeSolver::readControls()
{
    dict().readIfPresent("minIter", minIter_);
    dict().readIfPresent("maxIter", maxIter_);
    dict().readIfPresent("tolerance", tolerance_);
    dict().readIfPresent("relTol", relTolerance_);
}


Foam::scalar Foam::coupledIterativeSolver::normFactor
(
    const FieldField<Field, scalar>& x,
    const FieldField<Field, scalar>& b,
    const FieldField<Field, scalar>& Ax,
    FieldField<Field, scalar>& tmpField,
    const direction cmpt
) const
{
    typedef FieldField<Field, scalar> scalarFieldField;

    scalarFieldField pA(x.size());

    forAll (x, rowI)
    {
        // Calculate reference value of x
        scalar xRef = gAverage(x[rowI]);

        pA.set(rowI, new scalarField(x[rowI].size(), xRef));
    }

    // Calculate A.xRefField
    matrix_.Amul(tmpField, pA, bouCoeffs_, interfaces_, cmpt);

    // Calculate the normalisation factor
    return gSum(mag(Ax - tmpField) + mag(b - tmpField)) + lduMatrix::small_;
}


Foam::scalar Foam::coupledIterativeSolver::normFactor
(
    const FieldField<Field, scalar>& x,
    const FieldField<Field, scalar>& b,
    const direction cmpt
) const
{
    typedef FieldField<Field, scalar> scalarFieldField;

    scalarFieldField wA(x.size());
    scalarFieldField tmpField(x.size());

    forAll (x, rowI)
    {
        wA.set(rowI, new scalarField(x[rowI].size(), 0));
        tmpField.set(rowI, new scalarField(x[rowI].size()));
    }

    // Calculate A.x
    matrix_.Amul(wA, x, bouCoeffs_, interfaces_, cmpt);

    return normFactor(x, b, wA, tmpField, cmpt);
}


bool Foam::coupledIterativeSolver::stop
(
    coupledSolverPerformance& solverPerf
) const
{
    if (solverPerf.nIterations() < minIter_)
    {
        return false;
    }

    if
    (
        solverPerf.nIterations() >= maxIter_
     || solverPerf.checkConvergence(tolerance_, relTolerance_)
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
