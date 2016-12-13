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

template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix()
:
    Matrix<SquareMatrix<Type>, Type>()
{}


template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix(const label n)
:
    Matrix<SquareMatrix<Type>, Type>(n, n)
{}


template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix
(
    const label n,
    const Type& t
)
:
    Matrix<SquareMatrix<Type>, Type>(n, n, t)
{}


template<class Type>
inline Foam::SquareMatrix<Type>::SquareMatrix(Istream& is)
:
    Matrix<SquareMatrix<Type>, Type>(is)
{}


template<class Type>
inline Foam::autoPtr<Foam::SquareMatrix<Type> >
Foam::SquareMatrix<Type>::clone() const
{
    return autoPtr<SquareMatrix<Type> >(new SquareMatrix<Type>(*this));
}


// ************************************************************************* //
