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

template<class Type>
inline const Foam::Field< Foam::Field<Type> >&
Foam::correlationFunction<Type>::tZeroBuffers() const
{
    return tZeroBuffers_;
}


template<class Type>
inline Foam::scalar Foam::correlationFunction<Type>::duration() const
{
    return duration_;
}


template<class Type>
inline Foam::scalar Foam::correlationFunction<Type>::sampleInterval() const
{
    return sampleInterval_;
}


template<class Type>
inline Foam::scalar Foam::correlationFunction<Type>::averagingInterval() const
{
    return averagingInterval_;
}


template<class Type>
inline Foam::label Foam::correlationFunction<Type>::sampleSteps() const
{
    return sampleSteps_;
}


template<class Type>
inline Foam::label Foam::correlationFunction<Type>::measurandFieldSize() const
{
    return tZeroBuffers_[0].size();
}


// ************************************************************************* //
