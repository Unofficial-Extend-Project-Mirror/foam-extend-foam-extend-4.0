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

#include "includeEntry.H"
#include "dictionary.H"
#include "IFstream.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::functionEntries::includeEntry::typeName
(
    Foam::functionEntries::includeEntry::typeName_()
);

// Don't lookup the debug switch here as the debug switch dictionary
// might include includeEntry
Foam::debug::debugSwitch
Foam::functionEntries::includeEntry::debug
(
    "includeEntry",
    0
);

namespace Foam
{
namespace functionEntries
{
    addToMemberFunctionSelectionTable
    (
        functionEntry,
        includeEntry,
        execute,
        dictionaryIstream
    );

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        includeEntry,
        execute,
        primitiveEntryIstream
    );
}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::functionEntries::includeEntry::includeFileName
(
    Istream& is
)
{
    fileName fName(is);
    fName.expand();

    if (fName.size() && fName[0] != '/' && fName[1] != ':')
    {
        fName = fileName(is.name()).path()/fName;
    }

    return fName;
}


bool Foam::functionEntries::includeEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    IFstream fileStream(includeFileName(is));

    if (fileStream)
    {
        parentDict.read(fileStream);
        return true;
    }
    else
    {
        FatalIOErrorIn
        (
            "functionEntries::includeEntry::includeEntry"
            "(dictionary& parentDict,Istream& is)",
            is
        )   << "Cannot open include file " << fileStream.name()
            << " while reading dictionary " << parentDict.name()
            << exit(FatalIOError);

        return false;
    }
}

bool Foam::functionEntries::includeEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    IFstream fileStream(includeFileName(is));

    if (fileStream)
    {
        entry.read(parentDict, fileStream);
        return true;
    }
    else
    {
        FatalIOErrorIn
        (
            "functionEntries::includeEntry::includeEntry"
            "(dictionary& parentDict, primitiveEntry& entry, Istream& is)",
            is
        )   << "Cannot open include file " << fileStream.name()
            << " while reading dictionary " << parentDict.name()
            << exit(FatalIOError);

        return false;
    }
}

// ************************************************************************* //
