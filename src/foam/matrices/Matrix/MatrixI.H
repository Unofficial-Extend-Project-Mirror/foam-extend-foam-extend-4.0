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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Form, class Type>
inline Foam::Matrix<Form, Type>::Matrix()
:
    v_(NULL),
    n_(0),
    m_(0)
{}


template<class Form, class Type>
inline Foam::autoPtr<Foam::Matrix<Form, Type> >
Foam::Matrix<Form, Type>::clone() const
{
    return autoPtr<Matrix<Form, Type> >(new Matrix<Form, Type>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Form, class Type>
inline const Foam::Matrix<Form, Type>& Foam::Matrix<Form, Type>::null()
{
    return zero;
}


//- Return the number of rows
template<class Form, class Type>
inline Foam::label Foam::Matrix<Form, Type>::n() const
{
    return n_;
}


template<class Form, class Type>
inline Foam::label Foam::Matrix<Form, Type>::m() const
{
    return m_;
}


template<class Form, class Type>
inline Foam::label Foam::Matrix<Form, Type>::size() const
{
    return n_*m_;
}


template<class Form, class Type>
inline void Foam::Matrix<Form, Type>::checki(const label i) const
{
    if (!n_)
    {
        FatalErrorIn("Matrix<Form, Type>::checki(const label)")
            << "attempt to access element from zero sized row"
            << abort(FatalError);
    }
    else if (i<0 || i>=n_)
    {
        FatalErrorIn("Matrix<Form, Type>::checki(const label)")
            << "index " << i << " out of range 0 ... " << n_-1
            << abort(FatalError);
    }
}


template<class Form, class Type>
inline void Foam::Matrix<Form, Type>::checkj(const label j) const
{
    if (!m_)
    {
        FatalErrorIn("Matrix<Form, Type>::checkj(const label)")
            << "attempt to access element from zero sized column"
            << abort(FatalError);
    }
    else if (j<0 || j>=m_)
    {
        FatalErrorIn("Matrix<Form, Type>::checkj(const label)")
            << "index " << j << " out of range 0 ... " << m_-1
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Form, class Type>
inline Type* Foam::Matrix<Form, Type>::operator[](const label i)
{
#   ifdef FULLDEBUG
    checki(i);
#   endif
    return v_[i];
}


template<class Form, class Type>
inline const Type* Foam::Matrix<Form, Type>::operator[](const label i) const
{
#   ifdef FULLDEBUG
    checki(i);
#   endif
    return v_[i];
}


// ************************************************************************* //
