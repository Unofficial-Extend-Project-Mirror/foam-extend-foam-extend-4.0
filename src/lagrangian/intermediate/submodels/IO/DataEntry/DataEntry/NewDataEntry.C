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

#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::DataEntry<Type> > Foam::DataEntry<Type>::New
(
    const word& entryName,
    const dictionary& dict
)
{
    Istream& is(dict.lookup(entryName));

    word DataEntryType(is);

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(DataEntryType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("DataEntry<Type>::New(Istream&)")
            << "Unknown DataEntry type " << DataEntryType << " for DataEntry "
            << entryName << ". Constructor not in hash table" << nl << nl
            << "    Valid DataEntry types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<DataEntry<Type> >(cstrIter()(entryName, is));
}


// ************************************************************************* //
