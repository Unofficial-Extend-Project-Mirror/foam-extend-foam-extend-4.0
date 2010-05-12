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

#ifndef HashTable_C
#define HashTable_C

#include "HashTable.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
HashTable<T, Key, Hash>::HashTable(const label size)
:
    tableSize_(size),
    table_(NULL),
    nElmts_(0),
    endIter_(*this, NULL, 0),
    endConstIter_(*this, NULL, 0)
{
    if (tableSize_)
    {
        table_ = new hashedEntry*[tableSize_];
        for (label i=0; i<tableSize_; i++)
        {
            table_[i] = 0;
        }
    }
}


template<class T, class Key, class Hash>
HashTable<T, Key, Hash>::HashTable(const HashTable<T, Key, Hash>& ht)
:
    HashTableName(),
    tableSize_(ht.tableSize_),
    table_(NULL),
    nElmts_(0),
    endIter_(*this, NULL, 0),
    endConstIter_(*this, NULL, 0)
{
    if (tableSize_)
    {
        table_ = new hashedEntry*[tableSize_];

        for (label i=0; i<tableSize_; i++)
        {
            table_[i] = 0;
        }

        for (const_iterator iter = ht.begin(); iter != ht.end(); ++iter)
        {
            insert(iter.key(), *iter);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
HashTable<T, Key, Hash>::~HashTable()
{
    if (table_)
    {
        clear();
        delete[] table_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
bool HashTable<T, Key, Hash>::found(const Key& key) const
{
    if (tableSize_)
    {
        label ii = Hash()(key, tableSize_);

        for (hashedEntry* n=table_[ii]; n; n=n->next_)
        {
            if (key == n->key_) return true;
        }
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Info<< "HashTable<T, Key, Hash>::found(const Key& key) : "
            << "Entry " << key << " not found in hash table\n";
    }
#   endif

    return false;
}


template<class T, class Key, class Hash>
typename HashTable<T, Key, Hash>::iterator HashTable<T, Key, Hash>::find
(
    const Key& key
)
{
    if (tableSize_)
    {
        label ii = Hash()(key, tableSize_);
        hashedEntry* prev = 0;

        for (hashedEntry* n=table_[ii]; n; n=n->next_)
        {
            if (key == n->key_) return iterator(*this, n, ii);
            prev = n;
        }
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Info<< "HashTable<T, Key, Hash>::find(const Key& key) : "
            << "Entry " << key << " not found in hash table\n";
    }
#   endif

    return end();
}


template<class T, class Key, class Hash>
typename HashTable<T, Key, Hash>::const_iterator HashTable<T, Key, Hash>::find
(
    const Key& key
) const
{
    if (tableSize_)
    {
        label ii = Hash()(key, tableSize_);
        hashedEntry* prev = 0;

        for (hashedEntry* n=table_[ii]; n; n=n->next_)
        {
            if (key == n->key_) return const_iterator(*this, n, ii);
            prev = n;
        }
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Info<< "HashTable<T, Key, Hash>::find(const Key& key) const : "
            << "Entry " << key << " not found in hash table\n";
    }
#   endif

    return end();
}


// Return the table of contents
template<class T, class Key, class Hash>
List<Key> HashTable<T, Key, Hash>::toc() const
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
bool HashTable<T, Key, Hash>::set
(
    const Key& key,
    const T& newEntry,
    const bool protect
)
{
    if (tableSize_ == 0)
    {
        resize(2);
    }

    label ii = Hash()(key, tableSize_);
    hashedEntry* existing = 0;
    hashedEntry* prev = 0;

    for (hashedEntry* curr = table_[ii]; curr; curr = curr->next_)
    {
        if (key == curr->key_)
        {
            existing = curr;
            break;
        }
        prev = curr;
    }

    // not found, insert it at the head
    if (!existing)
    {
        table_[ii] = new hashedEntry(key, table_[ii], newEntry);
        nElmts_++;

        if (double(nElmts_)/tableSize_ > 0.8)
        {
#           ifdef FULLDEBUG
            if (debug)
            {
                Info<< "HashTable<T, Key, Hash>::set"
                    "(const Key& key, T newEntry) : "
                    "Doubling table size\n";
            }
#           endif

            resize(2*tableSize_);
        }
    }
    else if (protect)
    {
        // found - but protected from overwriting
        // this corresponds to the STL 'insert' convention
#       ifdef FULLDEBUG
        if (debug)
        {
            Info<< "HashTable<T, Key, Hash>::set"
                "(const Key& key, T newEntry, false) : "
                "Cannot insert " << key << " already in hash table\n";
        }
#       endif
        return false;
    }
    else
    {
        // found - overwrite existing entry
        // this corresponds to the Perl convention
        hashedEntry* elemPtr = new hashedEntry(key, existing->next_, newEntry);

        // replace existing element - within list or insert at the head
        if (prev)
        {
            prev->next_ = elemPtr;
        }
        else
        {
            table_[ii] = elemPtr;
        }

        delete existing;
    }

    return true;
}


template<class T, class Key, class Hash>
bool HashTable<T, Key, Hash>::erase(const iterator& cit)
{
    if (cit.elmtPtr_)    // note: endIter_ also has 0 elmtPtr_
    {
        iterator& it = const_cast<iterator&>(cit);

        // Search element before elmtPtr_
        hashedEntry* prevElmtPtr = 0;

        for (hashedEntry* n=table_[it.hashIndex_]; n; n=n->next_)
        {
            if (n == it.elmtPtr_)
            {
                break;
            }
            prevElmtPtr = n;
        }

        if (prevElmtPtr)
        {
            // Have element before elmtPtr
            prevElmtPtr->next_ = it.elmtPtr_->next_;
            delete it.elmtPtr_;
            it.elmtPtr_ = prevElmtPtr;
        }
        else
        {
            // elmtPtr is first element on SLlist
            table_[it.hashIndex_] = it.elmtPtr_->next_;
            delete it.elmtPtr_;

            // Search back for previous non-zero table entry
            while (--it.hashIndex_ >= 0 && !table_[it.hashIndex_])
            {}

            if (it.hashIndex_ >= 0)
            {
                // In table entry search for last element
                it.elmtPtr_ = table_[it.hashIndex_];

                while (it.elmtPtr_ && it.elmtPtr_->next_)
                {
                    it.elmtPtr_ = it.elmtPtr_->next_;
                }
            }
            else
            {
                // No previous found. Mark with special value which is
                // - not end()
                // - handled by operator++
                it.elmtPtr_ = reinterpret_cast<hashedEntry*>(this);
                it.hashIndex_ = -1;
            }
        }

        nElmts_--;

#       ifdef FULLDEBUG
        if (debug)
        {
            Info<< "HashTable<T, Key, Hash>::erase(iterator&) : "
                << "hashedEntry " << it.elmtPtr_->key_ << " removed.\n";
        }
#       endif

        return true;
    }
    else
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            Info<< "HashTable<T, Key, Hash>::erase(iterator&) : "
                << "cannot remove hashedEntry from hash table\n";
        }
#       endif

        return false;
    }
}


template<class T, class Key, class Hash>
bool HashTable<T, Key, Hash>::erase(const Key& key)
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
void HashTable<T, Key, Hash>::resize(const label newSize)
{
    if (newSize == tableSize_)
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            Info<< "HashTable<T, Key, Hash>::resize(const label newSize) : "
                << "new table size == old table size\n";
        }
#       endif

        return;
    }

    HashTable<T, Key, Hash>* newTable = new HashTable<T, Key, Hash>(newSize);

    for (const_iterator iter = begin(); iter != end(); ++iter)
    {
        newTable->insert(iter.key(), *iter);
    }

    label oldTableSize = tableSize_;
    tableSize_ = newTable->tableSize_;
    newTable->tableSize_ = oldTableSize;

    hashedEntry** oldTable = table_;
    table_ = newTable->table_;
    newTable->table_ = oldTable;

    delete newTable;
}


template<class T, class Key, class Hash>
void HashTable<T, Key, Hash>::clear()
{
    if (nElmts_)
    {
        for (label i=0; i<tableSize_; i++)
        {
            if (table_[i])
            {
                hashedEntry* n = table_[i];
                while(hashedEntry* next = n->next_)
                {
                    delete n;
                    n = next;
                }
                delete n;
                table_[i] = 0;
            }
        }
        nElmts_ = 0;
    }
}


template<class T, class Key, class Hash>
void HashTable<T, Key, Hash>::transfer(HashTable<T, Key, Hash>& ht)
{
    clear();
    delete[] table_;

    tableSize_ = ht.tableSize_;
    ht.tableSize_ = 0;

    table_ = ht.table_;
    ht.table_ = NULL;

    nElmts_ = ht.nElmts_;
    ht.nElmts_ = 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
void HashTable<T, Key, Hash>::operator=(const HashTable<T, Key, Hash>& ht)
{
    // Check for assignment to self
    if (this == &ht)
    {
        FatalErrorIn
        (
            "HashTable<T, Key, Hash>::operator="
            "(const HashTable<T, Key, Hash>&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }

    clear();

    for (const_iterator iter = ht.begin(); iter != ht.end(); ++iter)
    {
        insert(iter.key(), *iter);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

#include "HashTableIO.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
