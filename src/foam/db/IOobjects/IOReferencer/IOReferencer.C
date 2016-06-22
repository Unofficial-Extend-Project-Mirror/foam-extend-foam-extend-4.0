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

#include "IOReferencer.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class Type>
Foam::IOReferencer<Type>::IOReferencer
(
    const IOobject& io
)
:
    regIOobject(io),
    typePtr_(NULL)
{
    if
    (
        io.readOpt() != IOobject::NO_READ
     || io.writeOpt() != IOobject::NO_WRITE
    )
    {
        FatalErrorIn("IOReferencer<Type>::IOReferencer")
            << "IOReferencer can only be NO_READ, NO_WRITE."
            << abort(FatalError);
    }
}


template<class Type>
Foam::IOReferencer<Type>::IOReferencer
(
    const IOobject& io,
    Type* ptr
)
:
    regIOobject(io),
    typePtr_(ptr)
{
    if
    (
        (io.readOpt() != IOobject::NO_READ)
     || (io.writeOpt() != IOobject::NO_WRITE)
    )
    {
        FatalErrorIn("IOReferencer<Type>::IOReferencer")
            << "IOReferencer can only be NO_READ, NO_WRITE."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class Type>
Foam::IOReferencer<Type>::~IOReferencer()
{
    if (typePtr_)
    {
        delete typePtr_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::IOReferencer<Type>::writeData(Ostream& os) const
{
    // do nothing
    return os.good();
}


template<class Type>
const Type& Foam::IOReferencer<Type>::operator()() const
{
    if (!typePtr_)
    {
        FatalErrorIn
        (
            "IOReferencer<Type>::operator()"
        )
            << "Attempting to derefence a null typePtr - use IOReferencer::set"
            << "first."
            << abort(FatalError);
    }

    return *typePtr_;
}


template<class Type>
Type& Foam::IOReferencer<Type>::operator()()
{
    if (!typePtr_)
    {
        FatalErrorIn
        (
            "IOReferencer<Type>::operator()"
        )
            << "Attempting to derefence a null typePtr - use IOReferencer::set"
            << "first."
            << abort(FatalError);
    }

    return *typePtr_;
}


template<class Type>
void Foam::IOReferencer<Type>::set
(
    Type* ptr
)
{
    if (typePtr_)
    {
        delete typePtr_;
    }

    typePtr_ = ptr;
}


// ************************************************************************* //
