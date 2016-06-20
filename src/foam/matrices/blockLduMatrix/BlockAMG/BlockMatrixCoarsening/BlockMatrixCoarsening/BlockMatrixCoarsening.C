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
    BlockMatrixCoarsening

\*---------------------------------------------------------------------------*/

#include "BlockMatrixCoarsening.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::BlockMatrixCoarsening<Type> >
Foam::BlockMatrixCoarsening<Type>::New
(
    const word& coarseningType,
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict,
    const label groupSize,
    const label nCoarseCells
)
{
    typename matrixConstructorTable::iterator constructorIter =
        matrixConstructorTablePtr_->find(coarseningType);

    if (constructorIter == matrixConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "autoPtr<BlockMatrixCoarsening<Type> > "
            "BlockMatrixCoarsening<Type>::New\n"
            "(\n"
            "    const word& coarseningType,\n"
            "    const lduMatrix& matrix,\n"
            "    const label groupSize\n"
            "    const label nCoarseCells\n"
            ")"
        )   << "Unknown AMG coarsening type. " << coarseningType
            << endl << endl
            << "Valid AMG coarsening types are :" << endl
            << matrixConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<BlockMatrixCoarsening<Type> >
    (
        constructorIter()
        (
            matrix,
            dict,
            groupSize,
            nCoarseCells
        )
    );
}


// ************************************************************************* //
