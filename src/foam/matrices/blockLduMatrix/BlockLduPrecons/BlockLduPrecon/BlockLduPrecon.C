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

#include "BlockLduPrecon.H"
#include "blockNoPrecons.H"

template<class Type>
class BlockNoPrecon;

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::BlockLduPrecon<Type> > Foam::BlockLduPrecon<Type>::New
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict,
    const word keyword
)
{
    word preconName;

    // Handle primitive or dictionary entry
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

    if (matrix.diagonal())
    {
        // No preconditioning for the diagonal matrix
        return autoPtr<BlockLduPrecon<Type> >
        (
            new BlockNoPrecon<Type>
            (
                matrix,
                controls
            )
        );
    }
    else
    {
        typename dictionaryConstructorTable::iterator constructorIter =
            dictionaryConstructorTablePtr_->find(preconName);

        if (constructorIter == dictionaryConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "autoPtr<BlockLduPrecon> BlockLduPrecon::New\n"
                "(\n"
                "    const BlockLduMatrix<Type>& matrix,\n"
                "    const dictionary& dict\n"
                ")",
                dict
            )   << "Unknown matrix preconditioner " << preconName
                << endl << endl
                << "Valid matrix preconditioners are :" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<BlockLduPrecon<Type> >
        (
            constructorIter()
            (
                matrix,
                controls
            )
        );
    }
}


// ************************************************************************* //
