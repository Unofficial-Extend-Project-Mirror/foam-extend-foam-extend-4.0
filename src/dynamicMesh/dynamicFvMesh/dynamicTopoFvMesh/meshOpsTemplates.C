/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    meshOps

Description
    Various utility functions that perform mesh-related operations.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "Pair.H"
#include "meshOps.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace meshOps
{
// Parallel non-blocking send for fixed lists
template <class Type, unsigned Size>
inline void pWrite
(
    const label toID,
    const FixedList<Type, Size>& data
)
{
    OPstream::write
    (
        Pstream::blocking,
        toID,
        reinterpret_cast<const char*>(&data[0]),
        Size*sizeof(Type)
    );
}


// Parallel non-blocking receive for fixed lists
template <class Type, unsigned Size>
inline void pRead
(
    const label fromID,
    FixedList<Type, Size>& data
)
{
    IPstream::read
    (
        Pstream::blocking,
        fromID,
        reinterpret_cast<char*>(data.begin()),
        Size*sizeof(Type)
    );
}


// Parallel non-blocking send for lists
template <class Type>
inline void pWrite
(
    const label toID,
    const UList<Type>& data
)
{
    OPstream::write
    (
        Pstream::nonBlocking,
        toID,
        reinterpret_cast<const char*>(&data[0]),
        data.size()*sizeof(Type)
    );
}


// Parallel non-blocking receive for lists
template <class Type>
inline void pRead
(
    const label fromID,
    UList<Type>& data
)
{
    IPstream::read
    (
        Pstream::nonBlocking,
        fromID,
        reinterpret_cast<char*>(&data[0]),
        data.size()*sizeof(Type)
    );
}


// Utility method to size-up the list to include an item
template <class Type>
inline void sizeUpList
(
    const Type& item,
    List<Type>& list
)
{
    list.setSize(list.size() + 1, item);
}


// Utility method to size-down the list to remove an item
template <class Type>
inline void sizeDownList
(
    const Type& item,
    List<Type>& list
)
{
    label index = -1;

    if ((index = findIndex(list, item)) > -1)
    {
        meshOps::removeIndex(index, list);
    }
    else
    {
        FatalErrorIn
        (
            "inline void meshOps::sizeDownList"
            "(const Type& item, List<Type>& list)"
        )
            << nl << "Item: " << item
            << " was not found in list. " << nl
            << " List: " << nl << list
            << abort(FatalError);
    }
}


// Remove an item at a particular index in the list
template <class Type>
inline void removeIndex
(
    const label index,
    List<Type>& list
)
{
    // Create a new list
    List<Type> newList(list.size() - 1);

    // Copy individual items
    label n = 0;

    forAll(list, itemI)
    {
        if (itemI == index)
        {
            continue;
        }

        newList[n++] = list[itemI];
    }

    // Overwrite
    list.transfer(newList);
}


} // End namespace meshOps


} // End namespace Foam

// ************************************************************************* //
