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
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Classes * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline label StaticHashTable<T, Key, Hash>::size() const
{
    return nElmts_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline T& StaticHashTable<T, Key, Hash>::operator[](const Key& key)
{
    iterator iter = find(key);

    if (iter == end())
    {
        FatalErrorIn("StaticHashTable<T, Key, Hash>::operator[](const Key&)")
            << key << " not found in table.  Valid entries are "
            << toc()
            << exit(FatalError);
    }

    return *iter;
}

template<class T, class Key, class Hash>
inline const T& StaticHashTable<T, Key, Hash>::operator[](const Key& key) const
{
    const_iterator iter = find(key);

    if (iter == end())
    {
        FatalErrorIn
        (
            "StaticHashTable<T, Key, Hash>::operator[](const Key&) const"
        )   << key << " not found in table.  Valid entries are "
            << toc()
            << exit(FatalError);
    }

    return *iter;
}


// * * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::
Iterator
(
    TableRef curStaticHashTable,
    label hashIndex,
    label elementIndex
)
:
    curStaticHashTable_(curStaticHashTable),
    hashIndex_(hashIndex),
    elementIndex_(elementIndex)
{}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::
Iterator(const iterator& iter)
:
    curStaticHashTable_(iter.curStaticHashTable_),
    hashIndex_(iter.hashIndex_),
    elementIndex_(iter.elementIndex_)
{}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline void StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::
operator=(const iterator& iter)
{
    this->hashIndex_ = iter.hashIndex_;
    this->elementIndex_ = iter.elementIndex_;
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline bool StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::
operator==(const iterator& iter) const
{
    return hashIndex_ == iter.hashIndex_ && elementIndex_ == iter.elementIndex_;
}

template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline bool StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::
operator==(const const_iterator& iter) const
{
    return hashIndex_ == iter.hashIndex_ && elementIndex_ == iter.elementIndex_;
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline bool StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::
operator!=(const iterator& iter) const
{
    return !operator==(iter);
}

template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline bool StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::
operator!=(const const_iterator& iter) const
{
    return !operator==(iter);
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline TRef StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::
operator*()
{
    return curStaticHashTable_.objects_[hashIndex_][elementIndex_];
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline TRef StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::
operator()()
{
    return operator*();
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline
typename StaticHashTable<T, Key, Hash>::template Iterator
<
    TRef,
    TableRef
>&
StaticHashTable<T, Key, Hash>::Iterator
<
    TRef,
    TableRef
>::operator++()
{
    const List<T>& localObjects = curStaticHashTable_.objects_[hashIndex_];

    if (elementIndex_ == localObjects.size()-1)
    {
        elementIndex_ = 0;

        // Find first non-zero entry
        for
        (
            hashIndex_++;
            hashIndex_ < curStaticHashTable_.objects_.size();
            hashIndex_++
        )
        {
            if (curStaticHashTable_.objects_[hashIndex_].size() > 0)
            {
                break;
            }
        }
    }
    else
    {
        elementIndex_++;
    }

    return *this;
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline
typename StaticHashTable<T, Key, Hash>::template Iterator
<
    TRef,
    TableRef
>
StaticHashTable<T, Key, Hash>::Iterator
<
    TRef,
    TableRef
>::operator++
(
    int
)
{
    iterator tmp = *this;
    ++*this;
    return tmp;
}


template<class T, class Key, class Hash>
template<class TRef, class TableRef>
inline
const Key& StaticHashTable<T, Key, Hash>::Iterator<TRef, TableRef>::
key()
{
    return curStaticHashTable_.keys_[hashIndex_][elementIndex_];
}


template<class T, class Key, class Hash>
inline typename StaticHashTable<T, Key, Hash>::iterator
StaticHashTable<T, Key, Hash>::begin()
{
    // Find first non-empty entry
    forAll(keys_, i)
    {
        if (keys_[i].size() > 0)
        {
            return iterator(*this, i, 0);
        }
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Info<< "StaticHashTable is empty\n";
    }
#   endif

    return StaticHashTable<T, Key, Hash>::endIter_;
}


template<class T, class Key, class Hash>
inline const typename StaticHashTable<T, Key, Hash>::iterator&
StaticHashTable<T, Key, Hash>::end()
{
    return StaticHashTable<T, Key, Hash>::endIter_;
}


template<class T, class Key, class Hash>
inline typename StaticHashTable<T, Key, Hash>::const_iterator
StaticHashTable<T, Key, Hash>::begin() const
{
    // Find first non-empty entry
    forAll(keys_, i)
    {
        if (keys_[i].size() > 0)
        {
            return const_iterator(*this, i, 0);
        }
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Info<< "StaticHashTable is empty\n";
    }
#   endif

    return StaticHashTable<T, Key, Hash>::endConstIter_;
}


template<class T, class Key, class Hash>
inline const typename StaticHashTable<T, Key, Hash>::const_iterator&
StaticHashTable<T, Key, Hash>::end() const
{
    return StaticHashTable<T, Key, Hash>::endConstIter_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
