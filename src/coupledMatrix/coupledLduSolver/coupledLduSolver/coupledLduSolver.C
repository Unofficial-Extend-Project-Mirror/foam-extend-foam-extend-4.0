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
    Top level solver class which selects the solver relevant to the particular
    matrix structure.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "coupledLduMatrix.H"
#include "coupledLduSolver.H"
#include "coupledDiagonalSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineRunTimeSelectionTable
(
    coupledLduSolver,
    symMatrix
);

defineRunTimeSelectionTable
(
    coupledLduSolver,
    asymMatrix
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<coupledLduSolver> coupledLduSolver::New
(
    const word& fieldName,
    const coupledLduMatrix& matrix,
    const PtrList<FieldField<Field, scalar> >& bouCoeffs,
    const PtrList<FieldField<Field, scalar> >& intCoeffs,
    const lduInterfaceFieldPtrsListList& interfaces,
    const dictionary& dict
)
{
    word solverName(dict.lookup("solver"));

    if (matrix.diagonal())
    {
        return autoPtr<coupledLduSolver>
        (
            new coupledDiagonalSolver
            (
                fieldName,
                matrix,
                bouCoeffs,
                intCoeffs,
                interfaces
            )
        );
    }
    else if (matrix.symmetric())
    {
        symMatrixConstructorTable::iterator constructorIter =
            symMatrixConstructorTablePtr_->find(solverName);

        if (constructorIter == symMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "autoPtr<coupledLduSolver> coupledLduSolver::New\n"
                "(\n"
                "    const word& fieldName,\n"
                "    const coupledLduMatrix& matrix,\n"
                "    const PtrList<FieldField<Field, scalar> >& bouCoeffs,\n"
                "    const PtrList<FieldField<Field, scalar> >& intCoeffs,\n"
                "    const lduInterfaceFieldPtrsListList& interfaces,\n"
                "    const dictionary& dict\n"
                ")",
                dict
            )   << "Unknown symmetric matrix solver " << solverName
                << endl << endl
                << "Valid symmetric matrix solvers are :" << endl
                << symMatrixConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<coupledLduSolver>
        (
            constructorIter()
            (
                fieldName,
                matrix,
                bouCoeffs,
                intCoeffs,
                interfaces,
                dict
            )
        );
    }
    else if (matrix.asymmetric())
    {
        asymMatrixConstructorTable::iterator constructorIter =
            asymMatrixConstructorTablePtr_->find(solverName);

        if (constructorIter == asymMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "autoPtr<coupledLduSolver> coupledLduSolver::New\n"
                "(\n"
                "    const word& fieldName,\n"
                "    const coupledLduMatrix& matrix,\n"
                "    const PtrList<FieldField<Field, scalar> >& bouCoeffs,\n"
                "    const PtrList<FieldField<Field, scalar> >& intCoeffs,\n"
                "    const lduInterfaceFieldPtrsListList& interfaces,\n"
                "    const dictionary& dict\n"
                ")",
                dict
            )   << "Unknown asymmetric matrix solver " << solverName
                << endl << endl
                << "Valid asymmetric matrix solvers are :" << endl
                << asymMatrixConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<coupledLduSolver>
        (
            constructorIter()
            (
                fieldName,
                matrix,
                bouCoeffs,
                intCoeffs,
                interfaces,
                dict
            )
        );
    }
    else
    {
        FatalErrorIn
        (
            "autoPtr<coupledLduSolver> coupledLduSolver::New\n"
            "(\n"
            "    const word& fieldName,\n"
            "    const coupledLduMatrix& matrix,\n"
            "    const direction cmpt,\n"
            "    const PtrList<FieldField<Field, scalar> >& bouCoeffs,\n"
            "    const PtrList<FieldField<Field, scalar> >& intCoeffs,\n"
            "    const lduInterfaceFieldPtrsListList& interfaces,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "cannot solve incomplete matrix, "
               "no diagonal or off-diagonal coefficient"
            << exit(FatalError);

        return autoPtr<coupledLduSolver>(NULL);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
