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

#include "error.H"
#include "pTraits.H"
#include "Swap.H"

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

template<class T>
const Foam::UList<T> Foam::UList<T>::zero;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::UList<T>::UList()
:
    v_(0),
    size_(0)
{}


template<class T>
inline Foam::UList<T>::UList(T* __restrict__ v, label size)
:
    v_(v),
    size_(size)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const Foam::UList<T>& Foam::UList<T>::null()
{
    return zero;
}


// Reset in slicing.  HJ. 19/Oct/2008
template<class T>
inline void Foam::UList<T>::reset(T* __restrict__ v, label size)
{
    v_ = v;
    size_ = size;
}


template<class T>
inline Foam::label Foam::UList<T>::fcIndex(const label i) const
{
    return (i == size()-1 ? 0 : i+1);
}


template<class T>
inline Foam::label Foam::UList<T>::rcIndex(const label i) const
{
    return (i ? i-1 : size()-1);
}


// Check start is within valid range (0 ... size-1).
template<class T>
inline void Foam::UList<T>::checkStart(const label start) const
{
    if (start<0 || (start && start>=size_))
    {
        FatalErrorIn("UList<T>::checkStart(const label)")
            << "start " << start << " out of range 0 ... " << max(size_-1, 0)
            << abort(FatalError);
    }
}


// Check size is within valid range (0 ... size).
template<class T>
inline void Foam::UList<T>::checkSize(const label size) const
{
    if (size<0 || size>size_)
    {
        FatalErrorIn("UList<T>::checkSize(const label)")
            << "size " << size << " out of range 0 ... " << size_
            << abort(FatalError);
    }
}


// Check index i is within valid range (0 ... size-1).
template<class T>
inline void Foam::UList<T>::checkIndex(const label i) const
{
    if (!size_)
    {
        FatalErrorIn("UList<T>::checkIndex(const label)")
            << "attempt to access element from zero sized list"
            << abort(FatalError);
    }
    else if (i<0 || i>=size_)
    {
        FatalErrorIn("UList<T>::checkIndex(const label)")
            << "index " << i << " out of range 0 ... " << size_-1
            << abort(FatalError);
    }
}


template<class T>
inline T& Foam::UList<T>::first()
{
    return this->operator[](0);
}


template<class T>
inline const T& Foam::UList<T>::first() const
{
    return this->operator[](0);
}


template<class T>
inline T& Foam::UList<T>::last()
{
    return this->operator[](this->size()-1);
}


template<class T>
inline const T& Foam::UList<T>::last() const
{
    return this->operator[](this->size()-1);
}


template<class T>
inline const T* Foam::UList<T>::cdata() const
{
    return v_;
}


template<class T>
inline T* Foam::UList<T>::data()
{
    return v_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// element access
template<class T>
inline T& Foam::UList<T>::operator[](const label i)
{
#   ifdef FULLDEBUG
    checkIndex(i);
#   endif
    return v_[i];
}


namespace Foam
{

    // Template specialization for bool
    template<>
    inline const bool& UList<bool>::operator[](const label i) const
    {
        // lazy evaluation - return false for out-of-range
        if (i < size_)
        {
            return v_[i];
        }
        else
        {
            return pTraits<bool>::zero;
        }
    }

} // end of namespace Foam


// const element access
template<class T>
inline const T& Foam::UList<T>::operator[](const label i) const
{
#   ifdef FULLDEBUG
    checkIndex(i);
#   endif
    return v_[i];
}


// Allow cast to a const List<T>&
template<class T>
inline Foam::UList<T>::operator const Foam::List<T>&() const
{
    return *reinterpret_cast<const List<T>*>(this);
}


// * * * * * * * * * * * * * * STL Member Functions  * * * * * * * * * * * * //

template<class T>
inline typename Foam::UList<T>::iterator
Foam::UList<T>::begin()
{
    return v_;
}

template<class T>
inline typename Foam::UList<T>::const_iterator
Foam::UList<T>::begin() const
{
    return v_;
}

template<class T>
inline typename Foam::UList<T>::const_iterator
Foam::UList<T>::cbegin() const
{
    return v_;
}

template<class T>
inline typename Foam::UList<T>::iterator
Foam::UList<T>::end()
{
    return &v_[size_];
}

template<class T>
inline typename Foam::UList<T>::const_iterator
Foam::UList<T>::end() const
{
    return &v_[size_];
}

template<class T>
inline typename Foam::UList<T>::const_iterator
Foam::UList<T>::cend() const
{
    return &v_[size_];
}

template<class T>
inline typename Foam::UList<T>::iterator
Foam::UList<T>::rbegin()
{
    return &v_[size_-1];
}

template<class T>
inline typename Foam::UList<T>::const_iterator
Foam::UList<T>::rbegin() const
{
    return &v_[size_-1];
}

template<class T>
inline typename Foam::UList<T>::const_iterator
Foam::UList<T>::crbegin() const
{
    return &v_[size_-1];
}

template<class T>
inline typename Foam::UList<T>::iterator
Foam::UList<T>::rend()
{
    return &v_[-1];
}

template<class T>
inline typename Foam::UList<T>::const_iterator
Foam::UList<T>::rend() const
{
    return &v_[-1];
}

template<class T>
inline typename Foam::UList<T>::const_iterator
Foam::UList<T>::crend() const
{
    return &v_[-1];
}

template<class T>
inline Foam::label Foam::UList<T>::size() const
{
    return size_;
}


template<class T>
inline Foam::label Foam::UList<T>::max_size() const
{
    return labelMax;
}


template<class T>
inline bool Foam::UList<T>::empty() const
{
    return !size_;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T>
inline void Foam::reverse(UList<T>& ul, const label n)
{
    for (int i=0; i<n/2; i++)
    {
        Swap(ul[i], ul[n-1-i]);
    }
}

template<class T>
inline void Foam::reverse(UList<T>& ul)
{
    reverse(ul, ul.size());
}


// ************************************************************************* //
