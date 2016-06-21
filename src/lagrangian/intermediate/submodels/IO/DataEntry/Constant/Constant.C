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

#include "Constant.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Constant<Type>::Constant(const word& entryName, Istream& is)
:
    DataEntry<Type>(entryName),
    value_(is)
{}


template<class Type>
Foam::Constant<Type>::Constant(const Constant<Type>& cnst)
:
    DataEntry<Type>(cnst),
    value_(cnst.value_)
{}


template<>
Foam::Constant<Foam::label>::Constant(const word& entryName, Istream& is)
:
    DataEntry<label>(entryName),
    value_(readLabel(is))
{}


template<>
Foam::Constant<Foam::scalar>::Constant(const word& entryName, Istream& is)
:
    DataEntry<scalar>(entryName),
    value_(readScalar(is))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Constant<Type>::~Constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Constant<Type>::value(const scalar x) const
{
    return value_;
}


template<class Type>
Type Foam::Constant<Type>::integrate(const scalar x1, const scalar x2) const
{
    return (x2 - x1)*value_;
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

#include "ConstantIO.C"


// ************************************************************************* //
