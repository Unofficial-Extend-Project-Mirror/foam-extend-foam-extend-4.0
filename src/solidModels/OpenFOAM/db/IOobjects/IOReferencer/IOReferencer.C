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

#include "IOReferencer.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class Type>
Foam::IOReferencer<Type>::IOReferencer
(
    const IOobject& io
)
:
    regIOobject(io)
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


template<class Type>
Foam::IOReferencer<Type>::IOReferencer
(
    const IOobject& io,
    const Type& reference
)
:
    regIOobject(io),
    typePtr_(& reference)
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
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::IOReferencer<Type>::writeData(Ostream& os) const
{
    // do nothing
    return os.good();
}


template<class Type>
const Type&
    Foam::IOReferencer<Type>::operator()() const
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
    return * typePtr_;
}


template<class Type>
void Foam::IOReferencer<Type>::set
(
    const Type& reference
)
{
    typePtr_ = &reference;
}

// ************************************************************************* //
