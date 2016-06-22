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

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

template<class T>
const Foam::Xfer<T> Foam::Xfer<T>::zero;


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template<class T>
inline const Foam::Xfer<T>& Foam::Xfer<T>::null()
{
    return zero;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::Xfer<T>::Xfer(T* p)
:
    ptr_(p ? p : new T)
{}


template<class T>
inline Foam::Xfer<T>::Xfer(T& t, bool allowTransfer)
:
    ptr_(new T)
{
    if (allowTransfer)
    {
        ptr_->transfer(t);
    }
    else
    {
        ptr_->operator=(t);
    }
}


template<class T>
inline Foam::Xfer<T>::Xfer(const T& t)
:
    ptr_(new T)
{
    ptr_->operator=(t);
}


template<class T>
inline Foam::Xfer<T>::Xfer(const Xfer<T>& t)
:
    ptr_(new T)
{
    ptr_->transfer(*(t.ptr_));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T>
inline Foam::Xfer<T>::~Xfer()
{
    delete ptr_;
    ptr_ = 0;
}


// * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::Xfer<T>::operator=(T& t)
{
    ptr_->transfer(t);
}


template<class T>
inline void Foam::Xfer<T>::operator=(const Xfer<T>& t)
{
    // silently ignore attempted copy to self
    if (this != &t)
    {
        ptr_->transfer(*(t.ptr_));
    }
}


template<class T>
inline T& Foam::Xfer<T>::operator()() const
{
    return *ptr_;
}


template<class T>
inline T* Foam::Xfer<T>::operator->() const
{
    return ptr_;
}


// * * * * * * * * * * * * *  Helper Functions * * * * * * * * * * * * * * * //


template<class T>
inline Foam::Xfer<T> Foam::xferCopy(const T& t)
{
    return Foam::Xfer<T>(t);
}


template<class T>
inline Foam::Xfer<T> Foam::xferMove(T& t)
{
    return Foam::Xfer<T>(t, true);
}


template<class T>
inline Foam::Xfer<T> Foam::xferTmp(Foam::tmp<T>& tt)
{
    return Foam::Xfer<T>(tt(), tt.isTmp());
}


template<class To, class From>
inline Foam::Xfer<To> Foam::xferCopyTo(const From& t)
{
    Foam::Xfer<To> xf;
    xf() = t;
    return xf;
}


template<class To, class From>
inline Foam::Xfer<To> Foam::xferMoveTo(From& t)
{
    Foam::Xfer<To> xf;
    xf().transfer(t);
    return xf;
}


// ************************************************************************* //
