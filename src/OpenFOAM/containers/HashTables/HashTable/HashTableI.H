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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Classes * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline HashTable<T, Key, Hash>::hashedEntry::hashedEntry
(
    const Key& key,
    hashedEntry* next,
    const T& newEntry
)
:
    key_(key),
    next_(next),
    obj_(newEntry)
{}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline label HashTable<T, Key, Hash>::size() const
{
    return nElmts_;
}


template<class T, class Key, class Hash>
inline bool HashTable<T, Key, Hash>::insert(const Key& key, const T& newEntry)
{
    return set(key, newEntry, true);
}


template<class T, class Key, class Hash>
inline bool HashTable<T, Key, Hash>::set(const Key& key, const T& newEntry)
{
    return set(key, newEntry, false);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline T& HashTable<T, Key, Hash>::operator[](const Key& key)
{
    iterator iter = find(key);

    if (iter == end())
    {
        FatalErrorIn("HashTable<T, Key, Hash>::operator[](const Key&)")
            << key << " not found in table.  Valid entries are "
            << toc()
            << exit(FatalError);
    }

    return *iter;
}

template<class T, class Key, class Hash>
inline const T& HashTable<T, Key, Hash>::operator[](const Key& key) const
{
    const_iterator iter = find(key);

    if (iter == end())
    {
        FatalErrorIn("HashTable<T, Key, Hash>::operator[](const Key&) const")
            << key << " not found in table.  Valid entries are "
            << toc()
            << exit(FatalError);
    }

    return *iter;
}


template<class T, class Key, class Hash>
inline T& HashTable<T, Key, Hash>::operator()(const Key& key)
{
    iterator iter = find(key);

    if (iter == end())
    {
        insert(key, T());
        return *find(key);
    }
    else
    {
        return *iter;
    }
}


// * * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline HashTable<T, Key, Hash>::iterator::iterator
(
    HashTable<T, Key, Hash>& curHashTable,
    hashedEntry* elmt,
    label hashIndex
)
:
    curHashTable_(curHashTable),
    elmtPtr_(elmt),
    hashIndex_(hashIndex)
{}


template<class T, class Key, class Hash>
inline void HashTable<T, Key, Hash>::iterator::operator=(const iterator& iter)
{
    elmtPtr_ = iter.elmtPtr_;
    hashIndex_ = iter.hashIndex_;
}


template<class T, class Key, class Hash>
inline bool HashTable<T, Key, Hash>::iterator::operator==
(
    const iterator& iter
) const
{
    return elmtPtr_ == iter.elmtPtr_;
}


template<class T, class Key, class Hash>
inline bool HashTable<T, Key, Hash>::iterator::operator!=
(
    const iterator& iter
) const
{
    return elmtPtr_ != iter.elmtPtr_;
}


template<class T, class Key, class Hash>
inline bool HashTable<T, Key, Hash>::iterator::operator==
(
    const const_iterator& iter
) const
{
    return elmtPtr_ == iter.elmtPtr_;
}


template<class T, class Key, class Hash>
inline bool HashTable<T, Key, Hash>::iterator::operator!=
(
    const const_iterator& iter
) const
{
    return elmtPtr_ != iter.elmtPtr_;
}


template<class T, class Key, class Hash>
inline T& HashTable<T, Key, Hash>::iterator::operator*()
{
    return elmtPtr_->obj_;
}


template<class T, class Key, class Hash>
inline T& HashTable<T, Key, Hash>::iterator::operator()()
{
    return operator*();
}


template<class T, class Key, class Hash>
inline
typename HashTable<T, Key, Hash>::iterator&
HashTable<T, Key, Hash>::iterator::operator++()
{
    // Check for special value from erase. (sets hashIndex to -1)
    if (hashIndex_ >= 0)
    {
        // Do we have additional elements on the singly linked list?
        if (elmtPtr_ && elmtPtr_->next_)
        {
            elmtPtr_ = elmtPtr_->next_;
            return *this;
        }
    }

    // Step to the next table entry
    while
    (
        ++hashIndex_ < curHashTable_.tableSize_
     && !(elmtPtr_ = curHashTable_.table_[hashIndex_])
    )
    {}

    if (hashIndex_ == curHashTable_.tableSize_)
    {
        // make end iterator
        hashIndex_ = 0;
        elmtPtr_ = 0;
    }
    return *this;
}


template<class T, class Key, class Hash>
inline typename HashTable<T, Key, Hash>::iterator
HashTable<T, Key, Hash>::iterator::operator++
(
    int
)
{
    iterator tmp = *this;
    ++*this;
    return tmp;
}


template<class T, class Key, class Hash>
inline
const Key& HashTable<T, Key, Hash>::iterator::key()
{
    return elmtPtr_->key_;
}


template<class T, class Key, class Hash>
inline typename HashTable<T, Key, Hash>::iterator
HashTable<T, Key, Hash>::begin()
{
    label i = 0;

    while (table_ && !table_[i] && ++i < tableSize_)
    {}

    if (i == tableSize_)
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            Info<< "HashTable is empty\n";
        }
#       endif

        return HashTable<T, Key, Hash>::endIter_;
    }
    else
    {
        return iterator(*this, table_[i], i);
    }
}


template<class T, class Key, class Hash>
inline const typename HashTable<T, Key, Hash>::iterator&
HashTable<T, Key, Hash>::end()
{
    return HashTable<T, Key, Hash>::endIter_;
}


// * * * * * * * * * * * * * * * STL const_iterator * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline HashTable<T, Key, Hash>::const_iterator::const_iterator
(
    const HashTable<T, Key, Hash>& curHashTable,
    const hashedEntry* elmt,
    label hashIndex
)
:
    curHashTable_(curHashTable),
    elmtPtr_(elmt),
    hashIndex_(hashIndex)
{}


template<class T, class Key, class Hash>
inline HashTable<T, Key, Hash>::const_iterator::const_iterator
(
    const iterator& iter
)
:
    curHashTable_(iter.curHashTable_),
    elmtPtr_(iter.elmtPtr_),
    hashIndex_(iter.hashIndex_)
{}


template<class T, class Key, class Hash>
inline void HashTable<T, Key, Hash>::const_iterator::operator=
(
    const const_iterator& iter
)
{
    elmtPtr_ = iter.elmtPtr_;
    hashIndex_ = iter.hashIndex_;
}


template<class T, class Key, class Hash>
inline bool HashTable<T, Key, Hash>::const_iterator::operator==
(
    const const_iterator& iter
) const
{
    return elmtPtr_ == iter.elmtPtr_;
}


template<class T, class Key, class Hash>
inline bool HashTable<T, Key, Hash>::const_iterator::operator!=
(
    const const_iterator& iter
) const
{
    return elmtPtr_ != iter.elmtPtr_;
}


template<class T, class Key, class Hash>
inline bool HashTable<T, Key, Hash>::const_iterator::operator==
(
    const iterator& iter
) const
{
    return elmtPtr_ == iter.elmtPtr_;
}


template<class T, class Key, class Hash>
inline bool HashTable<T, Key, Hash>::const_iterator::operator!=
(
    const iterator& iter
) const
{
    return elmtPtr_ != iter.elmtPtr_;
}


template<class T, class Key, class Hash>
inline const T& HashTable<T, Key, Hash>::const_iterator::operator*()
{
    return elmtPtr_->obj_;
}

#ifndef __CINT__
template<class T, class Key, class Hash>
inline const T& HashTable<T, Key, Hash>::const_iterator::operator()()
{
    return operator*();
}
#endif

template<class T, class Key, class Hash>
inline
typename HashTable<T, Key, Hash>::const_iterator&
HashTable<T, Key, Hash>::const_iterator::operator++()
{
    if
    (
        !(elmtPtr_ = elmtPtr_->next_)
     && ++hashIndex_ < curHashTable_.tableSize_
     && !(elmtPtr_ = curHashTable_.table_[hashIndex_])
    )
    {
        while
        (
            ++hashIndex_ < curHashTable_.tableSize_
         && !(elmtPtr_ = curHashTable_.table_[hashIndex_])
        )
        {}
    }

    return *this;
}


template<class T, class Key, class Hash>
inline typename HashTable<T, Key, Hash>::const_iterator
HashTable<T, Key, Hash>::const_iterator::operator++
(
    int
)
{
    const_iterator tmp = *this;
    ++*this;
    return tmp;
}


template<class T, class Key, class Hash>
inline
const Key& HashTable<T, Key, Hash>::const_iterator::key()
{
    return elmtPtr_->key_;
}


template<class T, class Key, class Hash>
inline typename HashTable<T, Key, Hash>::const_iterator
HashTable<T, Key, Hash>::begin() const
{
    label i = 0;

    while (table_ && !table_[i] && ++i < tableSize_)
    {}

    if (i == tableSize_)
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            Info<< "HashTable is empty\n";
        }
#       endif

        return HashTable<T, Key, Hash>::endConstIter_;
    }
    else
    {
        return const_iterator(*this, table_[i], i);
    }
}


template<class T, class Key, class Hash>
inline const typename HashTable<T, Key, Hash>::const_iterator&
HashTable<T, Key, Hash>::end() const
{
    return HashTable<T, Key, Hash>::endConstIter_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
