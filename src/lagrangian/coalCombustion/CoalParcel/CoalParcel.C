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

#include "CoalParcel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::CoalParcel<ThermoType>::CoalParcel
(
    ReactingMultiphaseCloud<CoalParcel<ThermoType> >& owner,
    const vector& position,
    const label cellI
)
:
    ReactingMultiphaseParcel<CoalParcel<ThermoType> >(owner, position, cellI)
{}


template<class ThermoType>
Foam::CoalParcel<ThermoType>::CoalParcel
(
    ReactingMultiphaseCloud<CoalParcel<ThermoType> >& owner,
    const vector& position,
    const label cellI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const vector& U0,
    const scalarField& YMixture0,
    const scalarField& YGas0,
    const scalarField& YLiquid0,
    const scalarField& YSolid0,
    const typename
        ReactingMultiphaseParcel<CoalParcel<ThermoType> >::
        constantProperties& constProps
)
:
    ReactingMultiphaseParcel<CoalParcel<ThermoType> >
    (
        owner,
        position,
        cellI,
        typeId,
        nParticle0,
        d0,
        U0,
        YMixture0,
        YGas0,
        YLiquid0,
        YSolid0,
        constProps
    )
{}


template<class ThermoType>
Foam::CoalParcel<ThermoType>::CoalParcel
(
    const Cloud<CoalParcel<ThermoType> >& cloud,
    Istream& is,
    bool readFields
)
:
    ReactingMultiphaseParcel<CoalParcel<ThermoType> >(cloud, is, readFields)
{}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::CoalParcel<ThermoType>::~CoalParcel()
{}


// ************************************************************************* //
