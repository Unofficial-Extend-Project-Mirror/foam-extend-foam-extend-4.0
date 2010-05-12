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

Description

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct given size
template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const Foam::UList<T>& completeList,
    const Foam::List<label>& addresses
)
:
    completeList_(completeList),
    addresses_(addresses)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline Foam::label Foam::IndirectList<T>::size() const
{
    return addresses_.size();
}


template<class T>
inline const Foam::UList<T>& Foam::IndirectList<T>::
completeList() const
{
    return completeList_;
}


template<class T>
inline const Foam::List<Foam::label>& Foam::IndirectList<T>::addresses() const
{
    return addresses_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline const T& Foam::IndirectList<T>::operator[](const Foam::label i) const
{
    return completeList_[addresses_[i]];
}


// ************************************************************************* //
