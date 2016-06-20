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

#include "AverageIOField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::AverageIOField<Type>::AverageIOField
(
    const IOobject& io
)
:
    regIOobject(io)
{
    readStream(typeName) >> average_;
    readStream(typeName) >> static_cast<Field<Type>&>(*this);
    close();
}


template<class Type>
Foam::AverageIOField<Type>::AverageIOField
(
    const IOobject& io,
    const label size
)
:
    regIOobject(io),
    Field<Type>(size),
    average_(pTraits<Type>::zero)
{}


template<class Type>
Foam::AverageIOField<Type>::AverageIOField
(
    const IOobject& io,
    const Type& average,
    const Field<Type>& f
)
:
    regIOobject(io),
    Field<Type>(f),
    average_(average)
{
    if (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    {
        readStream(typeName)
            >> average_
            >> static_cast<Field<Type>&>(*this);
        close();
    }
}


template<class Type>
bool Foam::AverageIOField<Type>::writeData(Ostream& os) const
{
    os  << average_
        << token::NL
        << static_cast<const Field<Type>&>(*this);

    return os.good();
}


// ************************************************************************* //
