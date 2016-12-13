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

#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineRunTimeSelectionTable(lduPreconditioner, symMatrix);
    defineRunTimeSelectionTable(lduPreconditioner, asymMatrix);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::lduPreconditioner>
Foam::lduPreconditioner::New
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& dict,
    const word keyword
)
{
    word preconName;

    // handle primitive or dictionary entry
    const entry& e = dict.lookupEntry(keyword, false, false);
    if (e.isDict())
    {
        e.dict().lookup(keyword) >> preconName;
    }
    else
    {
        e.stream() >> preconName;
    }

    const dictionary& controls = e.isDict() ? e.dict() : dictionary::null;

    if (matrix.symmetric())
    {
        symMatrixConstructorTable::iterator constructorIter =
            symMatrixConstructorTablePtr_->find(preconName);

        if (constructorIter == symMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "lduPreconditioner::New\n"
                "(\n"
                "    const lduMatrix& matrix,\n"
                "    const FieldField<Field, scalar>& coupleBouCoeffs,\n"
                "    const FieldField<Field, scalar>& coupleIntCoeffs,\n"
                "    const lduInterfaceFieldPtrsList& interfaces,\n"
                "    const dictionary& dict,\n"
                "    const word keyword\n"
                ")",
                dict
            )   << "Unknown symmetric matrix preconditioner "
                << preconName << endl << endl
                << "Valid symmetric matrix preconditioners are :" << endl
                << symMatrixConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<lduPreconditioner>
        (
            constructorIter()
            (
                matrix,
                coupleBouCoeffs,
                coupleIntCoeffs,
                interfaces,
                controls
            )
        );
    }
    else if (matrix.asymmetric())
    {
        asymMatrixConstructorTable::iterator constructorIter =
            asymMatrixConstructorTablePtr_->find(preconName);

        if (constructorIter == asymMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "lduPreconditioner::New\n"
                "(\n"
                "    const lduMatrix& matrix,\n"
                "    const FieldField<Field, scalar>& coupleBouCoeffs,\n"
                "    const FieldField<Field, scalar>& coupleIntCoeffs,\n"
                "    const lduInterfaceFieldPtrsList& interfaces,\n"
                "    const dictionary& dict\n"
                ")",
                dict
            )   << "Unknown asymmetric matrix preconditioner "
                << preconName << endl << endl
                << "Valid asymmetric matrix preconditioners are :" << endl
                << asymMatrixConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<lduPreconditioner>
        (
            constructorIter()
            (
                matrix,
                coupleBouCoeffs,
                coupleIntCoeffs,
                interfaces,
                controls
            )
        );
    }
    else
    {
        FatalIOErrorIn
        (
            "lduPreconditioner::New\n"
            "(\n"
            "    const lduMatrix& matrix,\n"
            "    const FieldField<Field, scalar>& coupleBouCoeffs,\n"
            "    const FieldField<Field, scalar>& coupleIntCoeffs,\n"
            "    const lduInterfaceFieldPtrsList& interfaces,\n"
            "    const dictionary& dict\n"
            ")",
            dict
        )   << "cannot solve incomplete matrix, "
               "no diagonal or off-diagonal coefficient"
            << exit(FatalIOError);

        return autoPtr<lduPreconditioner>(NULL);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lduMatrix::preconditioner::preconditioner
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    matrix_(matrix),
    coupleBouCoeffs_(coupleBouCoeffs),
    coupleIntCoeffs_(coupleIntCoeffs),
    interfaces_(interfaces)
{}


// ************************************************************************* //
