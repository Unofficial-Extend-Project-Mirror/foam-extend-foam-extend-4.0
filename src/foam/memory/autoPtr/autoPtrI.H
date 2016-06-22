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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::autoPtr<T>::autoPtr(T* tPtr)
:
    ptr_(tPtr)
{}


template<class T>
inline Foam::autoPtr<T>::autoPtr(const autoPtr<T>& ap)
:
    ptr_(ap.ptr_)
{
    ap.ptr_ = 0;
}


template<class T>
inline Foam::autoPtr<T>::~autoPtr()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::autoPtr<T>::empty() const
{
    return !ptr_;
}


template<class T>
inline bool Foam::autoPtr<T>::valid() const
{
    return ptr_;
}


template<class T>
inline T* Foam::autoPtr<T>::ptr()
{
    T* ptr = ptr_;
    ptr_ = NULL;
    return ptr;
}


template<class T>
inline void Foam::autoPtr<T>::set(T* p)
{
    if (ptr_)
    {
        FatalErrorIn("void autoPtr<T>::set(T*)")
            << "object already allocated"
            << abort(FatalError);
    }

    ptr_ = p;
}


template<class T>
inline void Foam::autoPtr<T>::reset(T* p)
{
    if (ptr_)
    {
        delete ptr_;
    }

    ptr_ = p;
}


template<class T>
inline void Foam::autoPtr<T>::clear()
{
    reset(NULL);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline T& Foam::autoPtr<T>::operator()()
{
    if (!ptr_)
    {
        FatalErrorIn("T& autoPtr<T>::operator()()")
            << "object is not allocated"
            << abort(FatalError);
    }

    return *ptr_;
}


template<class T>
inline const T& Foam::autoPtr<T>::operator()() const
{
    if (!ptr_)
    {
        FatalErrorIn("const T& autoPtr<T>::operator()() const")
            << "object is not allocated"
            << abort(FatalError);
    }

    return *ptr_;
}


/*
template<class T>
inline T& Foam::autoPtr<T>::operator*()
{
    return operator()();
}


template<class T>
inline const T& Foam::autoPtr<T>::operator*() const
{
    return operator()();
}
*/


template<class T>
inline Foam::autoPtr<T>::operator const T&() const
{
    return operator()();
}


template<class T>
inline T* Foam::autoPtr<T>::operator->()
{
    if (!ptr_)
    {
        FatalErrorIn("autoPtr<T>::operator->()")
            << "object is not allocated"
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline const T* Foam::autoPtr<T>::operator->() const
{
    return const_cast<autoPtr<T>&>(*this).operator->();
}


template<class T>
inline void Foam::autoPtr<T>::operator=(const autoPtr<T>& ap)
{
    if (this != &ap)
    {
        reset(const_cast<autoPtr<T>&>(ap).ptr());
    }
}


// ************************************************************************* //
