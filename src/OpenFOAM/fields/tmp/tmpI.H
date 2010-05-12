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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline tmp<T>::tmp(T* tPtr)
:
    isTmp_(true),
    ptr_(tPtr),
    ref_(*tPtr)
{}


template<class T>
inline tmp<T>::tmp(const T& tRef)
:
    isTmp_(false),
    ptr_(NULL),
    ref_(tRef)
{}


template<class T>
inline tmp<T>::tmp(const tmp<T>& t)
:
    isTmp_(t.isTmp_),
    ptr_(t.ptr_),
    ref_(t.ref_)
{
    if (isTmp_)
    {
        if (ptr_)
        {
            ptr_->operator++();
        }
        else
        {
            FatalErrorIn("tmp<T>::tmp(const tmp<T>&)")
                << "attempted copy of a deallocated temporary"
                << abort(FatalError);
        }
    }
}


template<class T>
inline tmp<T>::~tmp()
{
    if (isTmp_ && ptr_)
    {
        if (ptr_->okToDelete())
        {
            delete ptr_;
            ptr_ = NULL;
        }
        else
        {
            ptr_->operator--();
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool tmp<T>::isTmp() const
{
    return isTmp_;
}


template<class T>
inline bool tmp<T>::valid() const
{
    return (!isTmp_ || (isTmp_ && ptr_));
}


template<class T>
inline T* tmp<T>::ptr() const
{
    if (isTmp_)
    {
         if (!ptr_)
         {
             FatalErrorIn("tmp<T>::ptr() const")
                 << "temporary deallocated"
                 << abort(FatalError);
         }

         T* ptr = ptr_;
         ptr_ = NULL;

         ptr->resetRefCount();

         return ptr;
    }
    else
    {
        return new T(ref_);
    }
}


template<class T>
inline void tmp<T>::clear()
{
    if (isTmp_ && ptr_) // && ptr_->okToDelete())
    {
        delete ptr_;
        ptr_ = NULL;
    }
}


template<class T>
inline void tmp<T>::clear() const
{
    const_cast<tmp<T>&>(*this).clear();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline T& tmp<T>::operator()()
{
    if (isTmp_)
    {
        if (!ptr_)
        {
            FatalErrorIn("T& tmp<T>::operator()()")
                << "temporary deallocated"
                << abort(FatalError);
        }

        return *ptr_;
    }
    else
    {
        // Note: const is cast away!
        // Perhaps there should be two refs, one for const and one for non const
        // and if the ref is actually const then you cannot return it here.
        // Another possibility would be to store a const ref and a flag to say
        // wether the tmp was constructed with a const or a non-const argument.
        return const_cast<T&>(ref_);
    }
}


template<class T>
inline const T& tmp<T>::operator()() const
{
    if (isTmp_)
    {
        if (!ptr_)
        {
            FatalErrorIn("const T& tmp<T>::operator()() const")
                << "temporary deallocated"
                << abort(FatalError);
        }

        return *ptr_;
    }
    else
    {
        return ref_;
    }
}


template<class T>
inline tmp<T>::operator const T&() const
{
    return operator()();
}


template<class T>
inline T* tmp<T>::operator->()
{
    if (isTmp_)
    {
         if (!ptr_)
         {
             FatalErrorIn("tmp<T>::operator->()")
                 << "temporary deallocated"
                 << abort(FatalError);
         }

         return ptr_;
    }
    else
    {
        return &const_cast<T&>(ref_);
    }
}


template<class T>
inline const T* tmp<T>::operator->() const
{
    return const_cast<tmp<T>&>(*this).operator->();
}


template<class T>
inline void tmp<T>::operator=(const tmp<T>& t)
{
    if (isTmp_ && ptr_)
    {
        if (ptr_->okToDelete())
        {
            delete ptr_;
            ptr_ = NULL;
        }
        else
        {
            ptr_->operator--();
        }
    }

    if (t.isTmp_)
    {
        isTmp_ = true;
        ptr_ = t.ptr_;

        if (ptr_)
        {
            ptr_->operator++();
        }
        else
        {
            FatalErrorIn("tmp<T>::operator=(const tmp<T>& t)")
                << "attempted copy of a deallocated temporary"
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn("tmp<T>::operator=(const tmp<T>& t)")
            << "attempted to assign to a const reference to constant object"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
