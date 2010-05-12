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

#include "ConeInjection.H"
#include "DataEntry.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::ConeInjection<CloudType>::nParcelsToInject
(
    const scalar time0,
    const scalar time1
) const
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return round((time1 - time0)*parcelsPerSecond_);
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::ConeInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
) const
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return volumeFlowRate_().integrate(time0, time1);
    }
    else
    {
        return 0.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConeInjection<CloudType>::ConeInjection
(
    const dictionary& dict,
    CloudType& owner
)
:
    InjectionModel<CloudType>(dict, owner, typeName),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    position_(this->coeffDict().lookup("position")),
    direction_(this->coeffDict().lookup("direction")),
    parcelsPerSecond_
    (
        readScalar(this->coeffDict().lookup("parcelsPerSecond"))
    ),
    volumeFlowRate_
    (
        DataEntry<scalar>::New
        (
            "volumeFlowRate",
            this->coeffDict()
        )
    ),
    Umag_
    (
        DataEntry<scalar>::New
        (
            "Umag",
            this->coeffDict()
        )
    ),
    thetaInner_
    (
        DataEntry<scalar>::New
        (
            "thetaInner",
            this->coeffDict()
        )
    ),
    thetaOuter_
    (
        DataEntry<scalar>::New
        (
            "thetaOuter",
            this->coeffDict()
        )
    ),
    parcelPDF_
    (
        pdf::New
        (
            this->coeffDict().subDict("parcelPDF"),
            owner.rndGen()
        )
    ),
    tanVec1_(vector::zero),
    tanVec2_(vector::zero)
{
    // Normalise direction vector
    direction_ /= mag(direction_);

    // Determine direction vectors tangential to direction
    vector tangent = vector::zero;
    scalar magTangent = 0.0;

    while (magTangent < SMALL)
    {
        vector v = this->owner().rndGen().vector01();

        tangent = v - (v & direction_)*direction_;
        magTangent = mag(tangent);
    }

    tanVec1_ = tangent/magTangent;
    tanVec2_ = direction_^tanVec1_;

    // Set total volume to inject
    this->volumeTotal_ = volumeFlowRate_().integrate(0.0, duration_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConeInjection<CloudType>::~ConeInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::ConeInjection<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::scalar Foam::ConeInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::vector Foam::ConeInjection<CloudType>::position
(
    const label,
    const scalar,
    const polyMeshInfo& meshInfo
)
{
    vector pos = position_;
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
                "Foam::vector Foam::ConeInjection<CloudType>::position"
            )   << "Could not determine 2-D case geometry" << nl
                << abort(FatalError);
        }
    }

    return pos;
}


template<class CloudType>
Foam::vector Foam::ConeInjection<CloudType>::velocity
(
    const label,
    const scalar time,
    const polyMeshInfo& meshInfo
)
{
    const scalar deg2Rad = mathematicalConstant::pi/180.0;

    scalar t = time - this->SOI_;
    scalar ti = thetaInner_().value(t);
    scalar to = thetaOuter_().value(t);
    scalar coneAngle = this->owner().rndGen().scalar01()*(to - ti) + ti;

    coneAngle *= deg2Rad;
    scalar alpha = sin(coneAngle);
    scalar dcorr = cos(coneAngle);
    scalar beta =
        2.0*mathematicalConstant::pi*this->owner().rndGen().scalar01();

    vector normal = alpha*(tanVec1_*cos(beta) + tanVec2_*sin(beta));
    vector dirVec = dcorr*direction_;
    dirVec += normal;

    // Remove empty component of velocity for slab cases
    if (meshInfo.caseIs2dSlab())
    {
        dirVec.component(meshInfo.emptyComponent()) = 0.0;
    }

    dirVec /= mag(dirVec);

    return Umag_().value(t)*dirVec;
}


template<class CloudType>
Foam::scalar Foam::ConeInjection<CloudType>::d0
(
    const label,
    const scalar
) const
{
    return parcelPDF_().sample();
}


// ************************************************************************* //
