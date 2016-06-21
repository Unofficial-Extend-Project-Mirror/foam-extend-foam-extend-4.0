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

#include "inputModeEntry.H"
#include "dictionary.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::functionEntries::inputModeEntry::typeName
(
    Foam::functionEntries::inputModeEntry::typeName_()
);

// Don't lookup the debug switch here as the debug switch dictionary
// might include inputModeEntries
Foam::debug::debugSwitch
Foam::functionEntries::inputModeEntry::debug
(
    "inputModeEntry",
    0
);

Foam::functionEntries::inputModeEntry::inputMode
    Foam::functionEntries::inputModeEntry::mode_(MERGE);

namespace Foam
{
namespace functionEntries
{
    addToMemberFunctionSelectionTable
    (
        functionEntry,
        inputModeEntry,
        execute,
        dictionaryIstream
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// we could combine this into execute() directly, but leave it here for now
void Foam::functionEntries::inputModeEntry::setMode(Istream& is)
{
    clear();

    word mode(is);
    if (mode == "merge" || mode == "default")
    {
        mode_ = MERGE;
    }
    else if (mode == "overwrite")
    {
        mode_ = OVERWRITE;
    }
    else if (mode == "protect")
    {
        mode_ = PROTECT;
    }
    else if (mode == "warn")
    {
        mode_ = WARN;
    }
    else if (mode == "error")
    {
        mode_ = FATALERROR;
    }
    else
    {
        WarningIn("Foam::functionEntries::inputModeEntry::setMode(Istream&)")
            << "unsupported input mode '" << mode
            << "' ... defaulting to 'merge'"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::inputModeEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    setMode(is);
    return true;
}


void Foam::functionEntries::inputModeEntry::clear()
{
    mode_ = MERGE;
}


bool Foam::functionEntries::inputModeEntry::merge()
{
    return mode_ == MERGE;
}


bool Foam::functionEntries::inputModeEntry::overwrite()
{
    return mode_ == OVERWRITE;
}


bool Foam::functionEntries::inputModeEntry::protect()
{
    return mode_ == PROTECT;
}

bool Foam::functionEntries::inputModeEntry::error()
{
    return mode_ == FATALERROR;
}


// ************************************************************************* //
