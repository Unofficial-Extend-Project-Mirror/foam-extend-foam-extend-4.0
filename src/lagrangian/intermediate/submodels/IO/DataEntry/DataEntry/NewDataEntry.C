/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

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
    word DataEntryType(dict.lookup(entryName));

    //    Info<< "Selecting DataEntry " << DataEntryType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(DataEntryType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "DataEntry<Type>::New(const dictionary&"
        )   << "Unknown DataEntry type "
            << DataEntryType << " for " << entryName
            << ", constructor not in hash table" << nl << nl
            << "    Valid DataEntry types are :" << nl
            << dictionaryConstructorTablePtr_->toc() << nl
            << exit(FatalError);
    }

    return autoPtr<DataEntry<Type> >(cstrIter()(entryName, dict));
}


// ************************************************************************* //
