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

#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
Foam::simpleMatrix<T>::simpleMatrix(const label mSize)
:
    scalarMatrix(mSize),
    source_(mSize, pTraits<T>::zero)
{}


template<class T>
Foam::simpleMatrix<T>::simpleMatrix
(
    const scalarMatrix& matrix,
    const Field<T>& source
)
:
    scalarMatrix(matrix),
    source_(source)
{}


template<class T>
Foam::simpleMatrix<T>::simpleMatrix(Istream& is)
:
    scalarMatrix(is),
    source_(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
Foam::Field<T> Foam::simpleMatrix<T>::solve() const
{
    // Since matrix and source are trashed during solution,
    // a copy is made.  HJ, 23/Dec/2008
    scalarMatrix tmpMatrix = *this;
    Field<T> sourceSol = source_;

    scalarMatrix::solve(tmpMatrix, sourceSol);

    return sourceSol;
}


template<class T>
Foam::Field<T> Foam::simpleMatrix<T>::LUsolve() const
{
    // Since matrix and source are trashed during solution,
    // a copy is made.  HJ, 23/Dec/2008
    scalarMatrix luMatrix = *this;
    Field<T> sourceSol = source_;

    scalarMatrix::LUsolve(luMatrix, sourceSol);

    return sourceSol;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::simpleMatrix<T>::operator=(const simpleMatrix<T>& m)
{
    if (this == &m)
    {
        FatalErrorIn("simpleMatrix<T>::operator=(const simpleMatrix<T>&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    if (n() != m.n())
    {
        FatalErrorIn("simpleMatrix<T>::operator=(const simpleMatrix<T>&)")
            << "Different size matrices"
            << abort(FatalError);
    }

    if (source_.size() != m.source_.size())
    {
        FatalErrorIn("simpleMatrix<T>::operator=(const simpleMatrix<T>&)")
            << "Different size source vectors"
            << abort(FatalError);
    }

    scalarMatrix::operator=(m);
    source_ = m.source_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class T>
Foam::simpleMatrix<T> Foam::operator+
(
    const simpleMatrix<T>& m1,
    const simpleMatrix<T>& m2
)
{
    return simpleMatrix<T>
    (
        static_cast<const scalarMatrix&>(m1)
      + static_cast<const scalarMatrix&>(m2),
        m1.source_ + m2.source_
    );
}


template<class T>
Foam::simpleMatrix<T> Foam::operator-
(
    const simpleMatrix<T>& m1,
    const simpleMatrix<T>& m2
)
{
    return simpleMatrix<T>
    (
        static_cast<const scalarMatrix&>(m1)
      - static_cast<const scalarMatrix&>(m2),
        m1.source_ - m2.source_
    );
}


template<class T>
Foam::simpleMatrix<T> Foam::operator*(const scalar s, const simpleMatrix<T>& m)
{
    return simpleMatrix<T>(s*m.matrix_, s*m.source_);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T>
Foam::Ostream& Foam::operator<<(Ostream& os, const simpleMatrix<T>& m)
{
    os << static_cast<const scalarMatrix&>(m) << nl << m.source_;
    return os;
}


// ************************************************************************* //
