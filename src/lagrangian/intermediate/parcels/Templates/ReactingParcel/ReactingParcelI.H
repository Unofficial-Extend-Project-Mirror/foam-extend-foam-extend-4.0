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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::ReactingParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& dict
)
:
    ThermoParcel<ParcelType>::constantProperties(dict),
    Tvap_(dimensionedScalar(dict.lookup("Tvap")).value()),
    Tbp_(dimensionedScalar(dict.lookup("Tbp")).value())
{}


template<class ParcelType>
inline Foam::ReactingParcel<ParcelType>::trackData::trackData
(
    ReactingCloud<ParcelType>& cloud,
    const constantProperties& constProps,
    const interpolation<scalar>& rhoInterp,
    const interpolation<vector>& UInterp,
    const interpolation<scalar>& muInterp,
    const interpolation<scalar>& TInterp,
    const interpolation<scalar>& CpInterp,
    const interpolation<scalar>& pInterp,
    const vector& g
)
:
    ThermoParcel<ParcelType>::trackData
    (
        cloud,
        constProps,
        rhoInterp,
        UInterp,
        muInterp,
        TInterp,
        CpInterp,
        g
    ),
    cloud_(cloud),
    constProps_(constProps),
    pInterp_(pInterp)
{}


template<class ParcelType>
inline Foam::ReactingParcel<ParcelType>::ReactingParcel
(
    ReactingCloud<ParcelType>& owner,
    const label typeId,
    const vector& position,
    const label celli,
    const scalar d0,
    const vector& U0,
    const scalar nParticle0,
    const scalarField& YGas0,
    const scalarField& YLiquid0,
    const scalarField& YSolid0,
    const scalarField& YMixture0,
    const constantProperties& constProps
)
:
    ThermoParcel<ParcelType>
    (
        owner,
        typeId,
        position,
        celli,
        d0,
        U0,
        nParticle0,
        constProps
    ),
    mass0_(0.0),
    YMixture_(YMixture0),
    YGas_(YGas0),
    YLiquid_(YLiquid0),
    YSolid_(YSolid0),
    pc_(0.0)
{
    // Set initial parcel mass
    mass0_ = this->mass();
}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::ReactingParcel<ParcelType>::constantProperties::Tvap() const
{
    return Tvap_;
}


template<class ParcelType>
inline Foam::scalar
Foam::ReactingParcel<ParcelType>::constantProperties::Tbp() const
{
    return Tbp_;
}


// * * * * * * * * * * * trackData Member Functions  * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::ReactingCloud<ParcelType>&
Foam::ReactingParcel<ParcelType>::trackData::cloud()
{
    return cloud_;
}


template<class ParcelType>
inline const typename Foam::ReactingParcel<ParcelType>::constantProperties&
Foam::ReactingParcel<ParcelType>::trackData::constProps() const
{
    return constProps_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ReactingParcel<ParcelType>::trackData::pInterp() const
{
    return pInterp_;
}


// * * * * * * * * * * ThermoParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::scalarField& Foam::ReactingParcel<ParcelType>::YGas() const
{
    return YGas_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::ReactingParcel<ParcelType>::YGas()
{
    return YGas_;
}


template<class ParcelType>
inline const Foam::scalarField& Foam::ReactingParcel<ParcelType>::YLiquid() const
{
    return YLiquid_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::ReactingParcel<ParcelType>::YLiquid()
{
    return YLiquid_;
}


template<class ParcelType>
inline const Foam::scalarField& Foam::ReactingParcel<ParcelType>::YSolid() const
{
    return YSolid_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::ReactingParcel<ParcelType>::YSolid()
{
    return YSolid_;
}


template<class ParcelType>
inline const Foam::scalarField& Foam::ReactingParcel<ParcelType>::YMixture() const
{
    return YMixture_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::ReactingParcel<ParcelType>::YMixture()
{
    return YMixture_;
}


template<class ParcelType>
inline Foam::scalar Foam::ReactingParcel<ParcelType>::mass0() const
{
    return mass0_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ReactingParcel<ParcelType>::mass0()
{
    return mass0_;
}


// ************************************************************************* //
