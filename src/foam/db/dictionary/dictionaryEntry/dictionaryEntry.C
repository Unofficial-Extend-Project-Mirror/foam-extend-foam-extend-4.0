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

#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dictionaryEntry::dictionaryEntry
(
    const keyType& key,
    const dictionary& parentDict,
    const dictionary& dict
)
:
    entry(key),
    dictionary(parentDict, dict)
{}


Foam::dictionaryEntry::dictionaryEntry
(
    const dictionary& parentDict,
    const dictionaryEntry& dictEnt
)
:
    entry(dictEnt),
    dictionary(parentDict, dictEnt)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::dictionaryEntry::startLineNumber() const
{
    if (size())
    {
        return first()->startLineNumber();
    }
    else
    {
        return -1;
    }
}

Foam::label Foam::dictionaryEntry::endLineNumber() const
{
    if (size())
    {
        return last()->endLineNumber();
    }
    else
    {
        return -1;
    }
}


Foam::ITstream& Foam::dictionaryEntry::stream() const
{
    FatalIOErrorIn("ITstream& dictionaryEntry::stream() const", *this)
        << "Attempt to return dictionary entry as a primitive"
        << abort(FatalIOError);

    return lookup("");
}


const Foam::dictionary& Foam::dictionaryEntry::dict() const
{
    return *this;
}


Foam::dictionary& Foam::dictionaryEntry::dict()
{
    return *this;
}


// ************************************************************************* //
