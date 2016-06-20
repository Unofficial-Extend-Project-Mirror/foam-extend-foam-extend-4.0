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

#include "dlLibraryTable.H"
#include "dictionary.H"
#include "fileNameList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TablePtr>
bool Foam::dlLibraryTable::open
(
    const dictionary& dict,
    const word& libsEntry,
    const TablePtr& tablePtr
)
{
    if (dict.found(libsEntry))
    {
        fileNameList libNames(dict.lookup(libsEntry));

        bool allOpened = (libNames.size() > 0);

        forAll(libNames, i)
        {
            const fileName& libName = libNames[i];

            label nEntries = 0;

            if (tablePtr)
            {
                nEntries = tablePtr->size();
            }

            bool opened = dlLibraryTable::open(libName);
            allOpened = opened && allOpened;

            if (opened && (!tablePtr || tablePtr->size() <= nEntries))
            {
                WarningIn
                (
                    "dlLibraryTable::open"
                    "(const dictionary& dict, const word& libsEntry, "
                    "const TablePtr tablePtr)"
                )   << "library " << libName
                    << " did not introduce any new entries"
                    << endl << endl;
            }
        }

        return allOpened;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
