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

#include "Matrix.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class T>
void Foam::Matrix<T>::allocate()
{
    if (n_ && m_)
    {
        v_ = new T*[n_];
        v_[0] = new T[n_*m_];

        for (register label i=1; i<n_; i++)
        {
            v_[i] = v_[i-1] + m_;
        }
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class T>
Foam::Matrix<T>::~Matrix()
{
    if (v_)
    {
        delete[] (v_[0]);
        delete[] v_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
const Foam::Matrix<T>& Foam::Matrix<T>::null()
{
    Matrix<T>* nullPtr = reinterpret_cast<Matrix<T>*>(NULL);
    return *nullPtr;
}


template<class T>
Foam::Matrix<T>::Matrix(const label n, const label m)
:
    v_(NULL),
    n_(n),
    m_(m)
{
    if (n_ < 0 || m_ < 0)
    {
        FatalErrorIn("Matrix<T>::Matrix(const label n, const label m)")
            << "bad n, m " << n_ << ", " << m_
            << abort(FatalError);
    }

    allocate();
}


template<class T>
Foam::Matrix<T>::Matrix(const label n, const label m, const T& a)
:
    v_(NULL),
    n_(n),
    m_(m)
{
    if (n_ < 0 || m_ < 0)
    {
        FatalErrorIn
        (
            "Matrix<T>::Matrix(const label n, const label m, const T&)"
        )   << "bad n, m " << n_ << ", " << m_
            << abort(FatalError);
    }

    allocate();

    if (v_)
    {
        T* v = v_[0];

        label nm = n_*m_;

        for (register label i=0; i<nm; i++)
        {
            v[i] = a;
        }
    }
}


template<class T>
Foam::Matrix<T>::Matrix(const Matrix<T>& a)
:
    v_(NULL),
    n_(a.n_),
    m_(a.m_)
{
    if (a.v_)
    {
        allocate();
        T* v = v_[0];
        const T* av = a.v_[0];

        label nm = n_*m_;
        for (register label i=0; i<nm; i++)
        {
            v[i] = av[i];
        }
    }
}

        
template<class T>
void Foam::Matrix<T>::clear()
{
    if (v_)
    {
        delete[] (v_[0]);
        delete[] v_;
    }
    n_ = 0;
    m_ = 0;
    v_ = NULL;
}


template<class T>
void Foam::Matrix<T>::transfer(Matrix<T>& a)
{
    clear();

    n_ = a.n_;
    a.n_ = 0;

    m_ = a.m_;
    a.m_ = 0;

    v_ = a.v_;
    a.v_ = NULL;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::Matrix<T>::operator=(const T& t)
{
    if (v_)
    {
        T* v = v_[0];

        label nm = n_*m_;
        for (register label i=0; i<nm; i++)
        {
            v[i] = t;
        }
    }
}


// Assignment operator. Takes linear time.
template<class T>
void Foam::Matrix<T>::operator=(const Matrix<T>& a)
{
    if (this == &a)
    {
        FatalErrorIn("Matrix<T>::operator=(const Matrix<T>&)")
            << "attempted assignment to self"
            << abort(FatalError);
    }

    if (n_ != a.n_ || m_ != a.m_)
    {
        clear();
        n_ = a.n_;
        m_ = a.m_;
        allocate();
    } 

    if (v_)
    {
        T* v = v_[0];
        const T* av = a.v_[0];

        label nm = n_*m_;
        for (register label i=0; i<nm; i++)
        {
            v[i] = av[i];
        }
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T>
const T& Foam::max(const Matrix<T>& a)
{
    label nm = a.n_*a.m_;

    if (nm)
    {
        label curMaxI = 0;
        const T* v = a.v_[0];

        for (register label i=1; i<nm; i++)
        {
            if (v[i] > v[curMaxI])
            {
                curMaxI = i;
            }
        }

        return v[curMaxI];
    }
    else
    {
        FatalErrorIn("max(const Matrix<T>&)")
            << "matrix is empty"
            << abort(FatalError);

        // Return in error to keep compiler happy
        return a[0][0];
    }
}


template<class T>
const T& Foam::min(const Matrix<T>& a)
{
    label nm = a.n_*a.m_;

    if (nm)
    {
        label curMinI = 0;
        const T* v = a.v_[0];

        for (register label i=1; i<nm; i++)
        {
            if (v[i] < v[curMinI])
            {
                curMinI = i;
            }
        }

        return v[curMinI];
    }
    else
    {
        FatalErrorIn("min(const Matrix<T>&)")
            << "matrix is empty"
            << abort(FatalError);

        // Return in error to keep compiler happy
        return a[0][0];
    }
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class T>
Foam::Matrix<T> Foam::operator-(const Matrix<T>& a)
{
    Matrix<T> na(a.n_, a.m_);

    if (a.n_ && a.m_)
    {
        T* nav = na.v_[0];
        const T* av = a.v_[0];

        label nm = a.n_*a.m_;
        for (register label i=0; i<nm; i++)
        {
            nav[i] = -av[i];
        }
    }

    return na;
}


template<class T>
Foam::Matrix<T> Foam::operator+(const Matrix<T>& a, const Matrix<T>& b)
{
    if (a.n_ != b.n_)
    {
        FatalErrorIn("Matrix<T>::operator+(const Matrix<T>&, const Matrix<T>&)")
            << "attempted add matrices with different number of rows: "
            << a.n_ << ", " << b.n_
            << abort(FatalError);
    }

    if (a.m_ != b.m_)
    {
        FatalErrorIn("Matrix<T>::operator+(const Matrix<T>&, const Matrix<T>&)")
            << "attempted add matrices with different number of columns: "
            << a.m_ << ", " << b.m_
            << abort(FatalError);
    }

    Matrix<T> ab(a.n_, a.m_);

    T* abv = ab.v_[0];
    const T* av = a.v_[0];
    const T* bv = b.v_[0];

    label nm = a.n_*a.m_;;
    for (register label i=0; i<nm; i++)
    {
        abv[i] = av[i] + bv[i];
    }

    return ab;
}


template<class T>
Foam::Matrix<T> Foam::operator-(const Matrix<T>& a, const Matrix<T>& b)
{
    if (a.n_ != b.n_)
    {
        FatalErrorIn("Matrix<T>::operator-(const Matrix<T>&, const Matrix<T>&)")
            << "attempted add matrices with different number of rows: "
            << a.n_ << ", " << b.n_
            << abort(FatalError);
    }

    if (a.m_ != b.m_)
    {
        FatalErrorIn("Matrix<T>::operator-(const Matrix<T>&, const Matrix<T>&)")
            << "attempted add matrices with different number of columns: "
            << a.m_ << ", " << b.m_
            << abort(FatalError);
    }

    Matrix<T> ab(a.n_, a.m_);

    T* abv = ab.v_[0];
    const T* av = a.v_[0];
    const T* bv = b.v_[0];

    label nm = a.n_*a.m_;;
    for (register label i=0; i<nm; i++)
    {
        abv[i] = av[i] - bv[i];
    }

    return ab;
}


template<class T>
Foam::Matrix<T> Foam::operator*(const scalar s, const Matrix<T>& a)
{
    Matrix<T> sa(a.n_, a.m_);

    if (a.n_ && a.m_)
    {
        T* sav = sa.v_[0];
        const T* av = a.v_[0];

        label nm = a.n_*a.m_;;
        for (register label i=0; i<nm; i++)
        {
            sav[i] = s*av[i];
        }
    }

    return sa;
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "MatrixIO.C"

// ************************************************************************* //
