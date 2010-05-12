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

#include "InjectionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InjectionModel<CloudType>::InjectionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:   dict_(dict),
    owner_(owner),
    coeffDict_(dict.subDict(type + "Coeffs")),
    SOI_(readScalar(coeffDict_.lookup("SOI"))),
    volumeTotal_(0.0),
    timeStep0_(0.0),
    nParcels_(0),
    volume_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::InjectionModel<CloudType>::~InjectionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const CloudType& Foam::InjectionModel<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
CloudType& Foam::InjectionModel<CloudType>::owner()
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary& Foam::InjectionModel<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary& Foam::InjectionModel<CloudType>::coeffDict() const
{
    return coeffDict_;
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::timeStart() const
{
    return SOI_;
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::volumeTotal() const
{
    return volumeTotal_;
}


template<class CloudType>
Foam::label Foam::InjectionModel<CloudType>::nParcels() const
{
    return nParcels_;
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::volume() const
{
    return volume_;
}


template<class CloudType>
Foam::scalar Foam::InjectionModel<CloudType>::volumeFraction() const
{
    return volume_/volumeTotal_;
}


template<class CloudType>
void Foam::InjectionModel<CloudType>::prepareForNextTimeStep
(
    const scalar time0,
    const scalar time1
)
{
    // Initialise values
    nParcels_ = 0;
    volume_ = 0.0;

    // Return if not started injection event
    if (time1 < SOI_)
    {
        timeStep0_ = time1;
        return;
    }

    // Make times relative to SOI
    scalar t0 = timeStep0_ - SOI_;
    scalar t1 = time1 - SOI_;

    // Number of parcels to inject
    nParcels_ = nParcelsToInject(t0, t1);

    // Volume of parcels to inject
    volume_ = volumeToInject(t0, t1);

    // Hold previous time if no parcels, but non-zero volume fraction
    if ((nParcels_ == 0) && (volume_ > 0.0))
    {
        // hold value of timeStep0_
    }
    else
    {
        // advance value of timeStep0_
        timeStep0_ = time1;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NewInjectionModel.C"

// ************************************************************************* //
