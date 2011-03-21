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

#include "error.H"

#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct given initial table size
template<class T, class Key, class Hash>
HashPtrTable<T, Key, Hash>::HashPtrTable(label size)
:
    HashTable<T*, Key, Hash>(size)
{}


// Construct as copy
template<class T, class Key, class Hash>
HashPtrTable<T, Key, Hash>::HashPtrTable(const HashPtrTable<T, Key, Hash>& ht)
:
    HashTable<T*, Key, Hash>()
{
    for (const_iterator iter = ht.begin(); iter != ht.end(); ++iter)
    {
        // Bug fix, Microsoft port.  HJ, 21/Mar/2011
        insert(iter.key(), iter()->clone().ptr());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
HashPtrTable<T, Key, Hash>::~HashPtrTable()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
T* HashPtrTable<T, Key, Hash>::remove(iterator& it)
{
    T* elemPtr = *it;
    HashTable<T*, Key, Hash>::erase(it);
    return elemPtr;
}


template<class T, class Key, class Hash>
bool HashPtrTable<T, Key, Hash>::erase(iterator& it)
{
    T* elemPtr = *it;

    if (HashTable<T*, Key, Hash>::erase(it))
    {
        if (elemPtr)
        {
            delete elemPtr;
        }

        return true;
    }
    else
    {
        return false;
    }
}


template<class T, class Key, class Hash>
void HashPtrTable<T, Key, Hash>::clear()
{
    for
    (
        iterator iter = this->begin();
        iter != this->end();
        ++iter
    )
    {
        delete *iter;
    }

    HashTable<T*, Key, Hash>::clear();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
void HashPtrTable<T, Key, Hash>::operator=
(
    const HashPtrTable<T, Key, Hash>& ht
)
{
    // Check for assignment to self
    if (this == &ht)
    {
        FatalErrorIn
        (
            "HashPtrTable<T, Key, Hash>::operator="
            "(const HashPtrTable<T, Key, Hash>&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }

    clear();

    for(const_iterator iter = ht.begin(); iter != ht.end(); ++iter)
    {
        // Bug fix, Microsoft port.  HJ, 21/Mar/2011
        insert(iter.key(), iter()->clone().ptr());
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

#include "HashPtrTableIO.C"

// ************************************************************************* //
