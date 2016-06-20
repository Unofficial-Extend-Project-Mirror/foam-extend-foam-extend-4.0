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

#include "NoDispersion.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoDispersion<CloudType>::NoDispersion
(
    const dictionary&,
    CloudType& owner
)
:
    DispersionModel<CloudType>(owner)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoDispersion<CloudType>::~NoDispersion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::NoDispersion<CloudType>::active() const
{
    return false;
}


template<class CloudType>
void Foam::NoDispersion<CloudType>::cacheFields(const bool)
{
// do nothing
}


template<class CloudType>
Foam::vector Foam::NoDispersion<CloudType>::update
(
    const scalar,
    const label,
    const vector&,
    const vector& Uc,
    vector&,
    scalar&
)
{
    // Do nothing
    return Uc;
}


// ************************************************************************* //
