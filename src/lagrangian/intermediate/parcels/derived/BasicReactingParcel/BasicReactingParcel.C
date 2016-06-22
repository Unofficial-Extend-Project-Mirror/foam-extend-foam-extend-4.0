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

#include "BasicReactingParcel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::BasicReactingParcel<ThermoType>::BasicReactingParcel
(
    ReactingCloud<BasicReactingParcel<ThermoType> >& owner,
    const vector& position,
    const label cellI
)
:
    ReactingParcel<BasicReactingParcel<ThermoType> >(owner, position, cellI)
{}


template<class ThermoType>
Foam::BasicReactingParcel<ThermoType>::BasicReactingParcel
(
    ReactingCloud<BasicReactingParcel<ThermoType> >& owner,
    const vector& position,
    const label cellI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const vector& U0,
    const scalarField& Y0,
    const typename ReactingParcel<BasicReactingParcel<ThermoType> >::
        constantProperties& constProps
)
:
    ReactingParcel<BasicReactingParcel<ThermoType> >
    (
        owner,
        position,
        cellI,
        typeId,
        nParticle0,
        d0,
        U0,
        Y0,
        constProps
    )
{}


template<class ThermoType>
Foam::BasicReactingParcel<ThermoType>::BasicReactingParcel
(
    const Cloud<BasicReactingParcel<ThermoType> >& cloud,
    Istream& is,
    bool readFields
)
:
    ReactingParcel<BasicReactingParcel<ThermoType> >(cloud, is, readFields)
{}


template<class ThermoType>
Foam::BasicReactingParcel<ThermoType>::BasicReactingParcel
(
    const BasicReactingParcel<ThermoType>& p
)
:
    ReactingParcel<BasicReactingParcel>(p)
{}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::BasicReactingParcel<ThermoType>::~BasicReactingParcel()
{}


// ************************************************************************* //
