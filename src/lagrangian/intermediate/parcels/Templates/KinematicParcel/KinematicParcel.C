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

#include "KinematicParcel.H"
#include "dimensionedConstants.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::updateCellQuantities
(
    TrackData& td,
    const scalar dt,
    const label celli
)
{
    rhoc_ = td.rhoInterp().interpolate(this->position(), celli);
    Uc_ = td.UInterp().interpolate(this->position(), celli);
    muc_ = td.muInterp().interpolate(this->position(), celli);

    // Apply dispersion components to carrier phase velocity
    Uc_ = td.cloud().dispersion().update
    (
        dt,
        celli,
        U_,
        Uc_,
        UTurb_,
        tTurb_
    );
}


template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::calcCoupled
(
    TrackData& td,
    const scalar dt,
    const label celli
)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Define local properties at beginning of timestep
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    const scalar mass0 = mass();
//    const vector U0 = U_;


    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // Initialise transfer terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = vector::zero;


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate velocity - update U
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar Cud = 0.0;
    const vector U1 = calcVelocity(td, dt, Cud, dUTrans);


    // ~~~~~~~~~~~~~~~~~~~~~~~
    // Accumulate source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~

    // Update momentum transfer
    td.cloud().UTrans()[celli] += nParticle_*dUTrans;

    // Accumulate coefficient to be applied in carrier phase momentum coupling
    td.cloud().UCoeff()[celli] += nParticle_*mass()*Cud;


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    this->U() = U1;
}


template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::calcUncoupled
(
    TrackData& td,
    const scalar dt,
    const label
)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // Initialise transfer terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = vector::zero;


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate velocity - update U
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar Cud = 0.0;
    this->U() = calcVelocity(td, dt, Cud, dUTrans);
}


template<class ParcelType>
template<class TrackData>
Foam::vector Foam::KinematicParcel<ParcelType>::calcVelocity
(
    TrackData& td,
    const scalar dt,
    scalar& Cud,
    vector& dUTrans
)
{
    // Correct carrier phase velocity for 2-D slab cases
    const polyMeshInfo& meshInfo = td.cloud().meshInfo();
    if (meshInfo.caseIs2dSlab())
    {
        Uc_.component(meshInfo.emptyComponent()) = 0.0;
    }

    // Return linearised term from drag model
    Cud = td.cloud().drag().Cu(U_ - Uc_, d_, rhoc_, rho_, muc_);


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Set new particle velocity
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Update velocity - treat as 3-D
    const scalar bp = 1.0/(Cud + VSMALL);
    const vector ap = Uc_/bp + rhoc_/rho_*td.g();

    vector Unew = td.cloud().UIntegrator().integrate(U_, dt, ap, bp);

//    Info<< "U_, Unew = " << U_ << ", " << Unew << endl;

    // Calculate the momentum transfer to the continuous phase
    dUTrans = -mass()*(Unew - U_);

    // Make corrections for 2-D cases
    if (meshInfo.caseIs2d())
    {
        if (meshInfo.caseIs2dSlab())
        {
            // Remove the slab normal parcel velocity component
            Unew.component(meshInfo.emptyComponent()) = 0.0;
            dUTrans.component(meshInfo.emptyComponent()) = 0.0;

            // Snap parcels to central plane
            this->position().component(meshInfo.emptyComponent()) =
                meshInfo.centrePoint().component(meshInfo.emptyComponent());
        }
        else if (meshInfo.caseIs2dWedge())
        {
            // Snap parcels to central plane
            this->position().component(meshInfo.emptyComponent()) = 0.0;
        }
        else
        {
            FatalErrorIn("void Foam::KinematicParcel::calcVelocity")
                << "Could not determine 2-D case geometry" << nl
                << abort(FatalError);
        }
    }

    return Unew;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
bool Foam::KinematicParcel<ParcelType>::move
(
    TrackData& td
)
{
    ParcelType& p = static_cast<ParcelType&>(*this);

    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = td.cloud().pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

    const scalar deltaT = mesh.time().deltaT().value();
    scalar tEnd = (1.0 - p.stepFraction())*deltaT;
    const scalar dtMax = tEnd;

    while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
    {
        // Set the Lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        // Remember which cell the Parcel is in
        // since this will change if a face is hit
        label celli = p.cell();

        dt *= p.trackToFace(p.position() + dt*U_, td);

        tEnd -= dt;
        p.stepFraction() = 1.0 - tEnd/deltaT;

        // Update cell based properties
        p.updateCellQuantities(td, dt, celli);

        if (td.cloud().coupled())
        {
            p.calcCoupled(td, dt, celli);
        }
        else
        {
            p.calcUncoupled(td, dt, celli);
        }

        if (p.onBoundary() && td.keepParticle)
        {
            if (p.face() > -1)
            {
                if
                (
                    isType<processorPolyPatch>
                        (pbMesh[p.patch(p.face())])
                )
                {
                    td.switchProcessor = true;
                }
            }
        }
    }

    return td.keepParticle;
}


template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::hitProcessorPatch
(
    const processorPolyPatch&,
    TrackData& td
)
{
    td.switchProcessor = true;
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::hitProcessorPatch
(
    const processorPolyPatch&,
    int&
)
{}


template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::hitWallPatch
(
    const wallPolyPatch& wpp,
    TrackData& td
)
{
    td.cloud().wallInteraction().correct(wpp, this->face(), U_);
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::hitWallPatch
(
    const wallPolyPatch&,
    int&
)
{}


template<class ParcelType>
template<class TrackData>
void Foam::KinematicParcel<ParcelType>::hitPatch
(
    const polyPatch&,
    TrackData& td
)
{
    td.keepParticle = false;
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::hitPatch
(
    const polyPatch&,
    int&
)
{}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::transformProperties
(
    const tensor& T
)
{
    Particle<ParcelType>::transformProperties(T);
    U_ = transform(T, U_);
}


template<class ParcelType>
void Foam::KinematicParcel<ParcelType>::transformProperties
(
    const vector& separation
)
{
    Particle<ParcelType>::transformProperties(separation);
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "KinematicParcelIO.C"

// ************************************************************************* //

