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

#include "List.H"
#include "ListLoopM.H"

#include "FixedList.H"
#include "PtrList.H"
#include "SLList.H"
#include "IndirectList.H"
#include "BiIndirectList.H"
#include "contiguous.H"
#include "UIndirectList.H"

#include <algorithm>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

// Construct with length specified
template<class T>
List<T>::List(const label s)
:
    UList<T>(NULL, s)
{
    if (this->size_ < 0)
    {
        FatalErrorIn("List<T>::List(const label size)")
            << "bad size " << this->size_
            << abort(FatalError);
    }

    if (this->size_)
    {
        this->v_ = new T[this->size_];
    }
    else
    {
        this->v_ = 0;
    }
}


// Construct with length and single value specified
template<class T>
List<T>::List(const label s, const T& a)
:
    UList<T>(NULL, s)
{
    if (this->size_ < 0)
    {
        FatalErrorIn("List<T>::List(const label size, const T a)")
            << "bad size " << this->size_
            << abort(FatalError);
    }

    if (this->size_)
    {
        this->v_ = new T[this->size_];

        List_ACCESS(T, (*this), vp);
        List_FOR_ALL((*this), i)
            List_ELEM((*this), vp, i) = a;
        List_END_FOR_ALL
    }
    else
    {
        this->v_ = 0;
    }
}


// Construct as copy
template<class T>
List<T>::List(const List<T>& a)
:
    UList<T>(NULL, a.size_)
{
    if (this->size_)
    {
        this->v_ = new T[this->size_];

#       ifdef USEMEMCPY
        if (contiguous<T>())
        {
            memcpy(this->v_, a.v_, this->byteSize());
        }
        else
#       endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            List_FOR_ALL((*this), i)
                List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
            List_END_FOR_ALL
        }
    }
    else
    {
        this->v_ = 0;
    }
}


// Construct as copy or re-use as specified.
template<class T>
List<T>::List(List<T>& a, bool reUse)
:
    UList<T>(NULL, a.size_)
{
    if (reUse)
    {
        this->v_ = a.v_;
        a.v_ = 0;
        a.size_ = 0;
    }
    else if (this->size_)
    {
        this->v_ = new T[this->size_];

#       ifdef USEMEMCPY
        if (contiguous<T>())
        {
            memcpy(this->v_, a.v_, this->byteSize());
        }
        else
#       endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            List_FOR_ALL((*this), i)
                List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
            List_END_FOR_ALL
        }
    }
    else
    {
        this->v_ = 0;
    }
}


// Construct given size and start and end iterators.
template<class T>
template<class InputIterator>
List<T>::List(InputIterator first, InputIterator last)
{
    label s = 0;
    for
    (
        InputIterator iter = first;
        iter != last;
        ++iter
    )
    {
        s++;
    }

    setSize(s);

    s = 0;

    for
    (
        InputIterator iter = first;
        iter != last;
        ++iter
    )
    {
        this->operator[](s++) = iter();
    }
}


// Construct as copy of FixedList<T, Size>
template<class T>
template<label Size>
List<T>::List(const FixedList<T, Size>& fl)
:
    UList<T>(NULL, Size)
{
    if (Size)
    {
        this->v_ = new T[this->size_];

        forAll(*this, i)
        {
            this->operator[](i) = fl[i];
        }
    }
    else
    {
        this->v_ = 0;
    }
}


// Construct as copy of PtrList<T>
template<class T>
List<T>::List(const PtrList<T>& sptrl)
:
    UList<T>(NULL, sptrl.size())
{
    if (this->size_)
    {
        this->v_ = new T[this->size_];

        forAll(*this, i)
        {
            this->operator[](i) = sptrl[i];
        }
    }
    else
    {
        this->v_ = 0;
    }
}


// Construct as copy of SLList<T>
template<class T>
List<T>::List(const SLList<T>& sll)
:
    UList<T>(NULL, sll.size())
{
    if (this->size_)
    {
        this->v_ = new T[this->size_];

        label i = 0;
        for
        (
            typename SLList<T>::const_iterator iter = sll.begin();
            iter != sll.end();
            ++iter
        )
        {
            this->operator[](i++) = iter();
        }
    }
    else
    {
        this->v_ = 0;
    }
}


// Construct as copy of IndirectList<T>
template<class T>
List<T>::List(const IndirectList<T>& idl)
:
    UList<T>(NULL, idl.size())
{
    if (this->size_)
    {
        this->v_ = new T[this->size_];

        forAll(*this, i)
        {
            this->operator[](i) = idl[i];
        }
    }
    else
    {
        this->v_ = 0;
    }
}

// Construct as copy of UIndirectList<T>
template<class T>
List<T>::List(const UIndirectList<T>& lst)
:
    UList<T>(NULL, lst.size())
{
    if (this->size_)
    {
        this->v_ = new T[this->size_];

        forAll(*this, i)
        {
            this->operator[](i) = lst[i];
        }
    }
}

// Construct as copy of BiIndirectList<T>
template<class T>
List<T>::List(const BiIndirectList<T>& idl)
:
    UList<T>(NULL, idl.size())
{
    if (this->size_)
    {
        this->v_ = new T[this->size_];

        forAll(*this, i)
        {
            this->operator[](i) = idl[i];
        }
    }
    else
    {
        this->v_ = 0;
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

// Destroy list elements
template<class T>
List<T>::~List()
{
    if (this->size_) delete[] this->v_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
const List<T>& List<T>::null()
{
    List<T>* nullPtr = reinterpret_cast<List<T>*>(NULL);
    return *nullPtr;
}


template<class T>
void List<T>::setSize(const label newSize)
{
    if (newSize < 0)
    {
        FatalErrorIn("List<T>::setSize(const label)")
            << "bad set size " << newSize
            << abort(FatalError);
    }

    if (newSize != this->size_)
    {
        if (newSize > 0)
        {
            T* nv = new T[label(newSize)];

            if (this->size_)
            {
                register label i = min(this->size_, newSize);

#               ifdef USEMEMCPY
                if (contiguous<T>())
                {
                    memcpy(nv, this->v_, i*sizeof(T));
                }
                else
#               endif
                {
                    register T* vv = &this->v_[i];
                    register T* av = &nv[i];
                    while (i--) *--av = *--vv;
                }

                delete[] this->v_;
            }

            this->size_ = newSize;
            this->v_ = nv;
        }
        else
        {
            clear();
        }
    }
}


template<class T>
void List<T>::setSize(const label newSize, const T& a)
{
    label oldSize = this->size_;
    this->setSize(newSize);

    if (newSize > oldSize)
    {
        register label i = newSize - oldSize;
        register T* vv = &this->v_[newSize];
        while (i--) *--vv = a;
    }
}


template<class T>
void List<T>::clear()
{
    if (this->size_) delete[] this->v_;
    this->size_ = 0;
    this->v_ = 0;
}


// Transfer the contents of the argument List into this List
// and anull the argument list
template<class T>
void List<T>::transfer(List<T>& a)
{
    if (this->size_) delete[] this->v_;

    this->size_ = a.size_;
    this->v_ = a.v_;

    a.size_ = 0;
    a.v_ = 0;
}


template<class T>
template<unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
void List<T>::transfer(DynamicList<T, SizeInc, SizeMult, SizeDiv>& a)
{
    // shrink the allocated space to the number of elements used
    a.shrink();

    if (this->size_) delete[] this->v_;

    this->size_ = a.size_;
    this->v_ = a.v_;

    a.size_ = 0;
    a.v_ = 0;
    a.nextFree_ = 0;
}


template<class T>
void sort(List<T>& a)
{
    std::sort(a.begin(), a.end());
}


template<class T, class Cmp>
void sort(List<T>& a, const Cmp& cmp)
{
    std::sort(a.begin(), a.end(), cmp);
}


template<class T>
void stableSort(List<T>& a)
{
    std::stable_sort(a.begin(), a.end());
}


template<class T, class Cmp>
void stableSort(List<T>& a, const Cmp& cmp)
{
    std::stable_sort(a.begin(), a.end(), cmp);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// Assignment to UList operator. Takes linear time.
template<class T>
void List<T>::operator=(const UList<T>& a)
{
    if (a.size_ != this->size_)
    {
        if (this->size_) delete[] this->v_;
        this->size_ = a.size_;
        if (this->size_) this->v_ = new T[this->size_];
    }

    if (this->size_)
    {
#       ifdef USEMEMCPY
        if (contiguous<T>())
        {
            memcpy(this->v_, a.v_, this->byteSize());
        }
        else
#       endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            List_FOR_ALL((*this), i)
                List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
            List_END_FOR_ALL
        }
    }
}


// Assignment operator. Takes linear time.
template<class T>
void List<T>::operator=(const List<T>& a)
{
    if (this == &a)
    {
        FatalErrorIn("List<T>::operator=(const List<T>&)")
            << "attempted assignment to self"
            << abort(FatalError);
    }

    operator=(static_cast<const UList<T>&>(a));
}


// Assignment operator. Takes linear time.
template<class T>
void List<T>::operator=(const SLList<T>& sll)
{
    if (sll.size() != this->size_)
    {
        if (this->size_) delete[] this->v_;
        this->size_ = sll.size();
        if (this->size_) this->v_ = new T[this->size_];
    }

    if (this->size_)
    {
        label i = 0;
        for
        (
            typename SLList<T>::const_iterator iter = sll.begin();
            iter != sll.end();
            ++iter
        )
        {
            this->operator[](i++) = iter();
        }
    }
}


// Assignment operator. Takes linear time.
template<class T>
void List<T>::operator=(const IndirectList<T>& idl)
{
    if (idl.size() != this->size_)
    {
        if (this->size_) delete[] this->v_;
        this->size_ = idl.size();
        if (this->size_) this->v_ = new T[this->size_];
    }

    if (this->size_)
    {
        forAll(*this, i)
        {
            this->operator[](i) = idl[i];
        }
    }
}

// Assignment operator. Takes linear time.
template<class T>
void List<T>::operator=(const UIndirectList<T>& lst)
{
    if (lst.size() != this->size_)
    {
        if (this->v_) delete[] this->v_;
        this->v_ = 0;
        this->size_ = lst.size();
        if (this->size_) this->v_ = new T[this->size_];
    }

    forAll(*this, i)
    {
        this->operator[](i) = lst[i];
    }
}

// Assignment operator. Takes linear time.
template<class T>
void List<T>::operator=(const BiIndirectList<T>& idl)
{
    if (idl.size() != this->size_)
    {
        if (this->size_) delete[] this->v_;
        this->size_ = idl.size();
        if (this->size_) this->v_ = new T[this->size_];
    }

    if (this->size_)
    {
        forAll(*this, i)
        {
            this->operator[](i) = idl[i];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "ListIO.C"

// ************************************************************************* //
