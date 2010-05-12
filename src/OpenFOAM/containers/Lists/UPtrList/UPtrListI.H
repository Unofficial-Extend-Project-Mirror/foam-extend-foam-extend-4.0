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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline label UPtrList<T>::size() const
{
    return ptrs_.size();
}


template<class T>
inline bool UPtrList<T>::set(const label i) const
{
    return ptrs_[i] != NULL;
}


template<class T>
inline T* UPtrList<T>::set(const label i, T* ptr)
{
    T* old = ptrs_[i];
    ptrs_[i] = ptr;
    return old;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
const T& UPtrList<T>::operator[](const label i) const
{
    if (!ptrs_[i])
    {
        FatalErrorIn("UPtrList::operator[] const")
            << "hanging pointer, cannot dereference"
            << abort(FatalError);
    }

    return *(ptrs_[i]);
}


template<class T>
T& UPtrList<T>::operator[](const label i)
{
    if (!ptrs_[i])
    {
        FatalErrorIn("UPtrList::operator[]")
            << "hanging pointer, cannot dereference"
            << abort(FatalError);
    }

    return *(ptrs_[i]);
}


template<class T>
const T* UPtrList<T>::operator()(const label i) const
{
    return ptrs_[i];
}


// * * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * //

template<class T>
inline UPtrList<T>::iterator::iterator(T** ptr)
:
    ptr_(ptr)
{}

template<class T>
inline bool UPtrList<T>::iterator::operator==(const iterator& iter) const
{
    return ptr_ == iter.ptr_;
}

template<class T>
inline bool UPtrList<T>::iterator::operator!=(const iterator& iter) const
{
    return ptr_ != iter.ptr_;
}

template<class T>
inline T& UPtrList<T>::iterator::operator*()
{
    return **ptr_;
}

template<class T>
inline T& UPtrList<T>::iterator::operator()()
{
    return operator*();
}

template<class T>
inline typename UPtrList<T>::iterator
UPtrList<T>::iterator::operator++()
{
    ++ptr_;
    return *this;
}

template<class T>
inline typename UPtrList<T>::iterator
UPtrList<T>::iterator::operator++(int)
{
    iterator tmp = *this;
    ++ptr_;
    return tmp;
}

template<class T>
inline typename UPtrList<T>::iterator
UPtrList<T>::iterator::operator--()
{
    --ptr_;
    return *this;
}

template<class T>
inline typename UPtrList<T>::iterator
UPtrList<T>::iterator::operator--(int)
{
    iterator tmp = *this;
    --ptr_;
    return tmp;
}

template<class T>
inline typename UPtrList<T>::iterator
UPtrList<T>::iterator::operator+=(label n)
{
    ptr_ += n;
    return *this;
}

template<class T>
inline typename UPtrList<T>::iterator
operator+(const typename UPtrList<T>::iterator& iter, label n)
{
    typename UPtrList<T>::iterator tmp = iter;
    return tmp += n;
}

template<class T>
inline typename UPtrList<T>::iterator
operator+(label n, const typename UPtrList<T>::iterator& iter)
{
    typename UPtrList<T>::iterator tmp = iter;
    return tmp += n;
}

template<class T>
inline typename UPtrList<T>::iterator
UPtrList<T>::iterator::operator-=(label n)
{
    ptr_ -= n;
    return *this;
}

template<class T>
inline typename UPtrList<T>::iterator
operator-(const typename UPtrList<T>::iterator& iter, label n)
{
    typename UPtrList<T>::iterator tmp = iter;
    return tmp -= n;
}

template<class T>
inline label operator-
(
    const typename UPtrList<T>::iterator& iter1,
    const typename UPtrList<T>::iterator& iter2
)
{
    return (iter1.ptr_ - iter2.ptr_)/sizeof(T*);
}

template<class T>
inline T& UPtrList<T>::iterator::operator[](label n)
{
    return *(*this + n);
}

template<class T>
inline bool UPtrList<T>::iterator::operator<(const iterator& iter) const
{
    return ptr_ < iter.ptr_;
}

template<class T>
inline bool UPtrList<T>::iterator::operator>(const iterator& iter) const
{
    return ptr_ > iter.ptr_;
}

template<class T>
inline bool UPtrList<T>::iterator::operator<=(const iterator& iter) const
{
    return ptr_ <= iter.ptr_;
}

template<class T>
inline bool UPtrList<T>::iterator::operator>=(const iterator& iter) const
{
    return ptr_ >= iter.ptr_;
}

template<class T>
inline typename UPtrList<T>::iterator UPtrList<T>::begin()
{
    return ptrs_.begin();
}

template<class T>
inline typename UPtrList<T>::iterator UPtrList<T>::end()
{
    return ptrs_.end();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
