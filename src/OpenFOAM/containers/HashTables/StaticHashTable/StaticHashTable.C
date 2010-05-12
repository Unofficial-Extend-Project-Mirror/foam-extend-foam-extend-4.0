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

#ifndef StaticHashTable_C
#define StaticHashTable_C

#include "StaticHashTable.H"
#include "List.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct given initial table size
template<class T, class Key, class Hash>
StaticHashTable<T, Key, Hash>::StaticHashTable(const label size)
:
    StaticHashTableName(),
    keys_(size),
    objects_(size),
    nElmts_(0),
    endIter_(*this, keys_.size(), 0),
    endConstIter_(*this, keys_.size(), 0)
{
    if (size < 1)
    {
        FatalErrorIn
        (
            "StaticHashTable<T, Key, Hash>::StaticHashTable(const label size)"
        )   << "Illegal size " << size << " for StaticHashTable."
            << " Minimum size is 1" << abort(FatalError);
    }
}


// Construct as copy
template<class T, class Key, class Hash>
StaticHashTable<T, Key, Hash>::StaticHashTable
(
    const StaticHashTable<T, Key, Hash>& ht
)
:
    StaticHashTableName(),
    keys_(ht.keys_),
    objects_(ht.objects_),
    nElmts_(ht.nElmts_),
    endIter_(*this, keys_.size(), 0),
    endConstIter_(*this, keys_.size(), 0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
StaticHashTable<T, Key, Hash>::~StaticHashTable()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
bool StaticHashTable<T, Key, Hash>::found(const Key& key) const
{
    label ii = Hash()(key, keys_.size());

    const List<Key>& localKeys = keys_[ii];

    forAll(localKeys, n)
    {
        if (localKeys[n] == key)
        {
            return true;
        }
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Pout<< "StaticHashTable<T, Key, Hash>::found(const Key& key) : "
            << "Entry " << key << " not found in hash table\n";
    }
#   endif

    return false;
}


template<class T, class Key, class Hash>
typename StaticHashTable<T, Key, Hash>::iterator
StaticHashTable<T, Key, Hash>::find
(
    const Key& key
)
{
    label ii = Hash()(key, keys_.size());

    const List<Key>& localKeys = keys_[ii];

    forAll(localKeys, n)
    {
        if (localKeys[n] == key)
        {
            return iterator(*this, ii, n);
        }
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Pout<< "StaticHashTable<T, Key, Hash>::find(const Key& key) : "
            << "Entry " << key << " not found in hash table\n";
    }
#   endif

    return end();
}


template<class T, class Key, class Hash>
typename StaticHashTable<T, Key, Hash>::const_iterator
StaticHashTable<T, Key, Hash>::find
(
    const Key& key
) const
{
    label ii = Hash()(key, keys_.size());

    const List<Key>& localKeys = keys_[ii];

    forAll(localKeys, n)
    {
        if (localKeys[n] == key)
        {
            return const_iterator(*this, ii, n);
        }
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Pout<< "StaticHashTable<T, Key, Hash>::find(const Key& key) const : "
            << "Entry " << key << " not found in hash table\n";
    }
#   endif

    return end();
}


// Return the table of contents
template<class T, class Key, class Hash>
List<Key> StaticHashTable<T, Key, Hash>::toc() const
{
    List<Key> tofc(nElmts_);

    label i = 0;

    for (const_iterator iter = begin(); iter != end(); ++iter)
    {
        tofc[i++] = iter.key();
    }

    return tofc;
}


template<class T, class Key, class Hash>
bool StaticHashTable<T, Key, Hash>::insert(const Key& key, const T& newEntry)
{
    label ii = Hash()(key, keys_.size());

    List<Key>& localKeys = keys_[ii];

    forAll(localKeys, n)
    {
        if (localKeys[n] == key)
        {
#           ifdef FULLDEBUG
            if (debug)
            {
                Pout<< "StaticHashTable<T, Key, Hash>::insert"
                       "(const Key& key, T newEntry) : "
                       "Cannot insert " << key << " already in hash table\n";
            }
#           endif

            return false;
        }
    }


    // Append.
    List<T>& localObjects = objects_[ii];

    label sz = localKeys.size();

    localKeys.setSize(sz+1);
    localObjects.setSize(sz+1);

    localKeys[sz] = key;
    localObjects[sz] = newEntry;

    nElmts_++;

    return true;
}


template<class T, class Key, class Hash>
bool StaticHashTable<T, Key, Hash>::erase(const iterator& it)
{
    if (it != end())
    {
        List<Key>& localKeys = keys_[it.hashIndex_];
        List<T>& localObjects = objects_[it.hashIndex_];

        // Copy down
        for (label i = it.elementIndex_+1; i < localKeys.size(); i++)
        {
            localKeys[i-1] = localKeys[i];
            localObjects[i-1] = localObjects[i];
        }
        localKeys.setSize(localKeys.size()-1);
        localObjects.setSize(localObjects.size()-1);

        nElmts_--;

#       ifdef FULLDEBUG
        if (debug)
        {
            Pout<< "StaticHashTable<T, Key, Hash>::erase(iterator&) : "
                << "hashedEntry removed.\n";
        }
#       endif

        return true;
    }
    else
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            Pout<< "StaticHashTable<T, Key, Hash>::erase(iterator&) : "
                << "cannot remove hashedEntry from hash table\n";
        }
#       endif

        return false;
    }
}


template<class T, class Key, class Hash>
bool StaticHashTable<T, Key, Hash>::erase(const Key& key)
{
    iterator it = find(key);

    if (it != end())
    {
        return erase(it);
    }
    else
    {
        return false;
    }
}


template<class T, class Key, class Hash>
void StaticHashTable<T, Key, Hash>::resize(const label newSize)
{
    if (newSize == keys_.size())
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            Pout<< "StaticHashTable<T, Key, Hash>::resize(const label) : "
                << "new table size == old table size\n";
        }
#       endif

        return;
    }

    if (newSize < 1)
    {
        FatalErrorIn
        (
            "StaticHashTable<T, Key, Hash>::StaticHashTable(const label size)"
        )   << "Illegal size " << newSize << " for StaticHashTable."
            << " Minimum size is 1" << abort(FatalError);
    }


    StaticHashTable<T, Key, Hash> newTable(newSize);

    for (iterator iter = begin(); iter != end(); ++iter)
    {
        newTable.insert(iter.key(), *iter);
    }

    transfer(newTable);

    // Adapt end() iterators
    endIter_.hashIndex_ = keys_.size();
    endConstIter_.hashIndex_ = keys_.size();
}


template<class T, class Key, class Hash>
void StaticHashTable<T, Key, Hash>::clear()
{
    forAll(keys_, ii)
    {
        keys_[ii].clear();
        objects_[ii].clear();
    }
}


template<class T, class Key, class Hash>
void StaticHashTable<T, Key, Hash>::transfer(StaticHashTable<T, Key, Hash>& ht)
{
    // Remove my existing elements
    clear();

    // Copy data from ht
    keys_.transfer(ht.keys_);
    objects_.transfer(ht.objects_);
    nElmts_ = ht.nElmts_;

    // Adapt end() iterators
    endIter_.hashIndex_ = keys_.size();
    endConstIter_.hashIndex_ = keys_.size();

    // Clear ht
    ht.nElmts_ = 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
void StaticHashTable<T, Key, Hash>::operator=
(
    const StaticHashTable<T, Key, Hash>& ht
)
{
    // Check for assignment to self
    if (this == &ht)
    {
        FatalErrorIn
        (
            "StaticHashTable<T, Key, Hash>::operator="
            "(const StaticHashTable<T, Key, Hash>&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }

    clear();

    for (const_iterator iter = ht.begin(); iter != ht.end(); ++iter)
    {
        insert(iter.key(), *iter);
    }

    // keys_.size() does not change so neither does end() iterator.
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

#include "StaticHashTableIO.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
