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

#include "UList.H"
#include "SLList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, Foam::label Size>
inline Foam::FixedList<T, Size>::FixedList()
{}

template<class T, Foam::label Size>
inline Foam::FixedList<T, Size>::FixedList(const T v[Size])
{
    for (register label i=0; i<Size; i++)
    {
        v_[i] = v[i];
    }
}

template<class T, Foam::label Size>
inline Foam::FixedList<T, Size>::FixedList(const T& t)
{
    for (register label i=0; i<Size; i++)
    {
        v_[i] = t;
    }
}

template<class T, Foam::label Size>
inline Foam::FixedList<T, Size>::FixedList(const UList<T>& ul)
{
    checkSize(ul.size());

    for (register label i=0; i<Size; i++)
    {
        v_[i] = ul[i];
    }
}

template<class T, Foam::label Size>
inline Foam::FixedList<T, Size>::FixedList(const SLList<T>& sll)
{
    checkSize(sll.size());

    register label i = 0;
    for
    (
        typename SLList<T>::const_iterator iter = sll.begin();
        iter != sll.end();
        ++iter
    )
    {
        operator[](i++) = iter();
    }
}

template<class T, Foam::label Size>
inline Foam::FixedList<T, Size>::FixedList(const FixedList<T, Size>& fl)
{
    for (register label i=0; i<Size; i++)
    {
        v_[i] = fl[i];
    }
}

template<class T, Foam::label Size>
inline Foam::autoPtr<Foam::FixedList<T, Size> >
Foam::FixedList<T, Size>::clone() const
{
    return autoPtr<FixedList<T, Size> >(new FixedList<T, Size>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, Foam::label Size>
inline Foam::label Foam::FixedList<T, Size>::fcIndex(const label i) const
{
    return (i == Size-1 ? 0 : i+1); 
}


template<class T, Foam::label Size>
inline Foam::label Foam::FixedList<T, Size>::rcIndex(const label i) const
{
    return (i == 0 ? Size-1 : i-1); 
}


// Check start is within valid range (0 ... size-1).
template<class T, Foam::label Size>
inline void Foam::FixedList<T, Size>::checkStart(const label start) const
{
    if (start<0 || (start && start>=Size))
    {
        FatalErrorIn("FixedList<T, Size>::checkStart(const label)")
            << "start " << start << " out of range 0 ... " << max(Size-1, 0)
            << abort(FatalError);
    }
}


// Check size is within valid range (0 ... size).
template<class T, Foam::label Size>
inline void Foam::FixedList<T, Size>::checkSize(const label size) const
{
    if (size<0 || size>Size)
    {
        FatalErrorIn("FixedList<T, Size>::checkSize(const label)")
            << "size " << size << " out of range 0 ... " << Size
            << abort(FatalError);
    }
}


// Check index i is within valid range (0 ... size-1).
template<class T, Foam::label Size>
inline void Foam::FixedList<T, Size>::checkIndex(const label i) const
{
    if (!Size)
    {
        FatalErrorIn("FixedList<T, Size>::checkIndex(const label)")
            << "attempt to access element from zero sized list"
            << abort(FatalError);
    }
    else if (i<0 || i>=Size)
    {
        FatalErrorIn("FixedList<T, Size>::checkIndex(const label)")
            << "index " << i << " out of range 0 ... " << Size-1
            << abort(FatalError);
    }
}


template<class T, Foam::label Size>
inline void Foam::FixedList<T, Size>::setSize(const label s)
{
#   ifdef FULLDEBUG
    checkSize(s);
#   endif
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// Return subscript-checked element access
template<class T, Foam::label Size>
inline T& Foam::FixedList<T, Size>::operator[](const label i)
{
#   ifdef FULLDEBUG
    checkIndex(i);
#   endif
    return v_[i];
}


// Return subscript-checked const element access
template<class T, Foam::label Size>
inline const T& Foam::FixedList<T, Size>::operator[](const label i) const
{
#   ifdef FULLDEBUG
    checkIndex(i);
#   endif
    return v_[i];
}


template<class T, Foam::label Size>
inline void Foam::FixedList<T, Size>::operator=(const T v[Size])
{
    for (register label i=0; i<Size; i++)
    {
        v_[i] = v[i];
    }
}

template<class T, Foam::label Size>
inline void Foam::FixedList<T, Size>::operator=(const UList<T>& ul)
{
    checkSize(ul.size());

    for (register label i=0; i<Size; i++)
    {
        v_[i] = ul[i];
    }
}

template<class T, Foam::label Size>
inline void Foam::FixedList<T, Size>::operator=(const SLList<T>& sll)
{
    checkSize(sll.size());

    register label i = 0;
    for
    (
        typename SLList<T>::const_iterator iter = sll.begin();
        iter != sll.end();
        ++iter
    )
    {
        operator[](i++) = iter();
    }
}

template<class T, Foam::label Size>
inline void Foam::FixedList<T, Size>::operator=(const T& t)
{
    for (register label i=0; i<Size; i++)
    {
        v_[i] = t;
    }
}


// * * * * * * * * * * * * * * STL Member Functions  * * * * * * * * * * * * //

template<class T, Foam::label Size>
inline typename Foam::FixedList<T, Size>::iterator
Foam::FixedList<T, Size>::begin()
{
    return v_;
}


template<class T, Foam::label Size>
inline typename Foam::FixedList<T, Size>::const_iterator
Foam::FixedList<T, Size>::begin() const
{
    return v_;
}


template<class T, Foam::label Size>
inline typename Foam::FixedList<T, Size>::iterator
Foam::FixedList<T, Size>::end()
{
    return &v_[Size];
}


template<class T, Foam::label Size>
inline typename Foam::FixedList<T, Size>::const_iterator
Foam::FixedList<T, Size>::end() const
{
    return &v_[Size];
}

template<class T, Foam::label Size>
inline typename Foam::FixedList<T, Size>::iterator
Foam::FixedList<T, Size>::rbegin()
{
    return &v_[Size-1];
}


template<class T, Foam::label Size>
inline typename Foam::FixedList<T, Size>::const_iterator
Foam::FixedList<T, Size>::rbegin() const
{
    return &v_[Size-1];
}


template<class T, Foam::label Size>
inline typename Foam::FixedList<T, Size>::iterator
Foam::FixedList<T, Size>::rend()
{
    return &v_[-1];
}


template<class T, Foam::label Size>
inline typename Foam::FixedList<T, Size>::const_iterator
Foam::FixedList<T, Size>::rend() const
{
    return &v_[-1];
}


template<class T, Foam::label Size>
inline Foam::label Foam::FixedList<T, Size>::size() const
{
    return Size;
}


template<class T, Foam::label Size>
inline Foam::label Foam::FixedList<T, Size>::max_size() const
{
    return Size;
}


template<class T, Foam::label Size>
inline bool Foam::FixedList<T, Size>::empty() const
{
    return false;
}


#ifndef __CINT__
template<class T, Foam::label Size>
template<class HashT>
inline Foam::FixedList<T, Size>::Hash<HashT>::Hash()
{}

//- Rotating Hash. From http://burtleburtle.net/bob/hash/doobs.html.
template<class T, Foam::label Size>
template<class HashT>
inline Foam::label Foam::FixedList<T, Size>::Hash<HashT>::operator()
(
    const FixedList<T, Size>& fl
) const
{
    static const label farbit(8*sizeof(label)-4);

    label val = Size;

    for (register int i=0; i<Size; i++)
    {
        val = (val<<4)^(val>>farbit)^HashT()(fl[i]);
    }

    return val;
}

template<class T, Foam::label Size>
template<class HashT>
inline Foam::label Foam::FixedList<T, Size>::Hash<HashT>::operator()
(
    const FixedList<T, Size>& fl,
    const label tableSize
) const
{
    return mag(operator()(fl)) % tableSize;
}
#endif


// ************************************************************************* //
