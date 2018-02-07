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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::regIOobject::ownedByRegistry() const
{
    return ownedByRegistry_;
}


inline void Foam::regIOobject::store()
{
    ownedByRegistry_ = true;
}


template<class Type>
inline Type& Foam::regIOobject::store(Type* tPtr)
{
    if (!tPtr)
    {
        FatalErrorIn("Type& regIOobject::store(Type*)")
            << "object deallocated"
            << abort(FatalError);
    }

    tPtr->regIOobject::ownedByRegistry_ = true;

    return *tPtr;
}


template<class Type>
inline Type& Foam::regIOobject::store(autoPtr<Type>& atPtr)
{
    Type* tPtr = atPtr.ptr();

    if (!tPtr)
    {
        FatalErrorIn
            (
                "Type& regIOobject::store(autoPtr<Type>&)"
            )   << "object deallocated"
                << abort(FatalError);
    }

    tPtr->regIOobject::ownedByRegistry_ = true;

    return *tPtr;
}


inline void Foam::regIOobject::release()
{
    ownedByRegistry_ = false;
}


inline void Foam::regIOobject::deregister()
{
    registered_ = false;
}


inline Foam::label Foam::regIOobject::eventNo() const
{
    return eventNo_;
}


inline Foam::label& Foam::regIOobject::eventNo()
{
    return eventNo_;
}


// ************************************************************************* //
