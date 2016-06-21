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
    BlockLduSmoother

Description
    Block LDU matrix smoother virtual base class

\*---------------------------------------------------------------------------*/

#include "BlockLduSmoother.H"

template<class Type>
class BlockNoSmoother;

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::BlockLduSmoother<Type> > Foam::BlockLduSmoother<Type>::New
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict,
    const word keyword
)
{
    word smootherName;

    // Handle primitive or dictionary entry
    const entry& e = dict.lookupEntry(keyword, false, false);
    if (e.isDict())
    {
        e.dict().lookup(keyword) >> smootherName;
    }
    else
    {
        e.stream() >> smootherName;
    }

    const dictionary& controls = e.isDict() ? e.dict() : dictionary::null;

    typename dictionaryConstructorTable::iterator constructorIter =
        dictionaryConstructorTablePtr_->find(smootherName);

    if (constructorIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "autoPtr<BlockLduSmoother> BlockLduSmoother::New\n"
            "(\n"
            "    const BlockLduMatrix<Type>& matrix,\n"
            "    const dictionary& dict,\n"
            "    const word keyword\n"
            ")",
            dict
        )   << "Unknown matrix smoother " << smootherName
            << endl << endl
            << "Valid matrix smoothers are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<BlockLduSmoother<Type> >
    (
        constructorIter()
        (
            matrix,
            controls
        )
    );
}


// ************************************************************************* //
