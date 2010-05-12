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

#include "CompositionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CompositionModel<CloudType>::CompositionModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:   dict_(dict),
    owner_(owner),
    coeffDict_(dict.subDict(type + "Coeffs")),
    carrierThermo_(owner.carrierThermo()),
    gases_(owner.gases()),
    liquids_
    (
        liquidMixture::New
        (
            owner.mesh().objectRegistry::lookupObject<dictionary>
            (
                "thermophysicalProperties"
            )
        )
    ),
    solids_
    (
        solidMixture::New
        (
            owner.mesh().objectRegistry::lookupObject<dictionary>
            (
                "thermophysicalProperties"
            )
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CompositionModel<CloudType>::~CompositionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const CloudType& Foam::CompositionModel<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary& Foam::CompositionModel<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary& Foam::CompositionModel<CloudType>::coeffDict() const
{
    return coeffDict_;
}


template<class CloudType>
const Foam::hCombustionThermo&
Foam::CompositionModel<CloudType>::carrierThermo() const
{
    return carrierThermo_;
}


template<class CloudType>
const Foam::PtrList<Foam::specieReactingProperties>&
Foam::CompositionModel<CloudType>::gases() const
{
    return gases_;
}


template<class CloudType>
const Foam::liquidMixture& Foam::CompositionModel<CloudType>::liquids() const
{
    return liquids_();
}


template<class CloudType>
const Foam::solidMixture& Foam::CompositionModel<CloudType>::solids() const
{
    return solids_();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NewCompositionModel.C"

// ************************************************************************* //

