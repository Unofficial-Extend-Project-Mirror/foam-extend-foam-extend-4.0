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

#include "Constant.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Constant<Type>::Constant
(
    const word& entryName,
    const dictionary& dict
)
:
    DataEntry<Type>(typeName, entryName, dict),
    value_(this->dict_.lookup("value"))
{}


template<>
Foam::Constant<Foam::label>::Constant
(
    const word& entryName,
    const dictionary& dict
)
:
    DataEntry<label>(typeName, entryName, dict),
    value_(readLabel(this->dict_.lookup("value")))
{}


template<>
Foam::Constant<Foam::scalar>::Constant
(
    const word& entryName,
    const dictionary& dict
)
:
    DataEntry<scalar>(typeName, entryName, dict),
    value_(readScalar(this->dict_.lookup("value")))
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


// ************************************************************************* //
