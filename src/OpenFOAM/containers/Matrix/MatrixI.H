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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::Matrix<T>::Matrix()
:
    v_(NULL),
    n_(0),
    m_(0)
{}


template<class T>
inline Foam::autoPtr<Foam::Matrix<T> > Foam::Matrix<T>::clone() const
{
    return autoPtr<Matrix<T> >(new Matrix<T>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return the number of rows
template<class T>
inline Foam::label Foam::Matrix<T>::n() const
{
    return n_;
}


//- Return the number of columns
template<class T>
inline Foam::label Foam::Matrix<T>::m() const
{
    return m_;
}


//- Return the number of columns
template<class T>
inline Foam::label Foam::Matrix<T>::size() const
{
    return n_*m_;
}


// Check index i is within valid range (0 ... n-1).
template<class T>
inline void Foam::Matrix<T>::checki(const label i) const
{
    if (!n_)
    {
        FatalErrorIn("Matrix<T>::checki(const label)")
            << "attempt to access element from zero sized row"
            << abort(FatalError);
    }
    else if (i<0 || i>=n_)
    {
        FatalErrorIn("Matrix<T>::checki(const label)")
            << "index " << i << " out of range 0 ... " << n_-1
            << abort(FatalError);
    }
}


// Check index j is within valid range (0 ... n-1).
template<class T>
inline void Foam::Matrix<T>::checkj(const label j) const
{
    if (!m_)
    {
        FatalErrorIn("Matrix<T>::checkj(const label)")
            << "attempt to access element from zero sized column"
            << abort(FatalError);
    }
    else if (j<0 || j>=m_)
    {
        FatalErrorIn("Matrix<T>::checkj(const label)")
            << "index " << j << " out of range 0 ... " << m_-1
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// Return subscript-checked element access
template<class T>
inline T* Foam::Matrix<T>::operator[](const label i)
{
#   ifdef FULLDEBUG
    checki(i);
#   endif
    return v_[i];
}


// Return subscript-checked const element access
template<class T>
inline const T* Foam::Matrix<T>::operator[](const label i) const
{
#   ifdef FULLDEBUG
    checki(i);
#   endif
    return v_[i];
}


// ************************************************************************* //
