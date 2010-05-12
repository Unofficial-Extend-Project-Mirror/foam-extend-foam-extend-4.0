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

#include "ManualInjection.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::ManualInjection<CloudType>::nParcelsToInject
(
    const scalar time0,
    const scalar time1
) const
{
    if ((0.0 >= time0) && (0.0 < time1))
    {
        return positions_.size();
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::ManualInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
) const
{
    // All parcels introduced at SOI
    if ((0.0 >= time0) && (0.0 < time1))
    {
        return this->volumeTotal_;
    }
    else
    {
        return 0.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ManualInjection<CloudType>::ManualInjection
(
    const dictionary& dict,
    CloudType& owner
)
:
    InjectionModel<CloudType>(dict, owner, typeName),
    positionsFile_(this->coeffDict().lookup("positionsFile")),
    positions_
    (
        IOobject
        (
            positionsFile_,
            owner.db().time().constant(),
            owner.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    diameters_(positions_.size()),
    U0_(this->coeffDict().lookup("U0")),
    parcelPDF_
    (
        pdf::New
        (
            this->coeffDict().subDict("parcelPDF"),
            owner.rndGen()
        )
    )
{
    // Construct parcel diameters
    forAll(diameters_, i)
    {
        diameters_[i] = parcelPDF_->sample();
    }

    // Determine volume of particles to inject
    this->volumeTotal_ = sum(pow(diameters_, 3))
        *mathematicalConstant::pi/6.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ManualInjection<CloudType>::~ManualInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::ManualInjection<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::scalar Foam::ManualInjection<CloudType>::timeEnd() const
{
    // Not used
    return 0.0;
}


template<class CloudType>
Foam::vector Foam::ManualInjection<CloudType>::position
(
    const label iParcel,
    const scalar time,
    const polyMeshInfo& meshInfo
)
{
    vector pos = positions_[iParcel];
    if (meshInfo.caseIs2d())
    {
        if (meshInfo.caseIs2dWedge())
        {
            pos.component(meshInfo.emptyComponent()) = 0.0;
        }
        else if (meshInfo.caseIs2dSlab())
        {
            pos.component(meshInfo.emptyComponent()) =
                meshInfo.centrePoint().component(meshInfo.emptyComponent());
        }
        else
        {
            FatalErrorIn
            (
                "Foam::vector Foam::ManualInjection<CloudType>::position"
            )   << "Could not determine 2-D case geometry" << nl
                << abort(FatalError);
        }
    }

    return pos;
}


template<class CloudType>
Foam::vector Foam::ManualInjection<CloudType>::velocity
(
    const label,
    const scalar,
    const polyMeshInfo& meshInfo
)
{
    vector vel = U0_;
    if (meshInfo.caseIs2dSlab())
    {
        vel.component(meshInfo.emptyComponent()) =
            meshInfo.centrePoint().component(meshInfo.emptyComponent());
    }

    return vel;
}


template<class CloudType>
Foam::scalar Foam::ManualInjection<CloudType>::d0
(
    const label iParcel,
    const scalar
) const
{
    return diameters_[iParcel];
}


// ************************************************************************* //
