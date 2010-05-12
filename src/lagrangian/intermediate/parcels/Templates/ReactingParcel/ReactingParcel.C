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

#include "ReactingParcel.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::updateCellQuantities
(
    TrackData& td,
    const scalar dt,
    const label celli
)
{
    ThermoParcel<ParcelType>::updateCellQuantities(td, dt, celli);

    pc_ = td.pInterp().interpolate(this->position(), celli);
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::calcCoupled
(
    TrackData& td,
    const scalar dt,
    const label celli
)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Define local properties at beginning of timestep
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const vector U0 = this->U_;
    const scalar mass0 = this->mass();
    const scalar np0 = this->nParticle_;
    const scalar T0 = this->T_;
    const scalar cp0 = this->cp_;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // Initialise transfer terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = vector::zero;

    // Enthalpy transfer from the particle to the carrier phase
    scalar dhTrans = 0.0;

    // Mass transfer from particle to carrier phase
    // - components exist in particle already
    scalarList dMassMT(td.cloud().gases().size(), 0.0);

    // Mass transfer due to surface reactions from particle to carrier phase
    // - components do not necessarily exist in particle already
    scalarList dMassSR(td.cloud().gases().size(), 0.0);

    // Total mass lost from particle due to surface reactions
    // - sub-model will adjust component mass fractions
    scalar dMassMTSR = 0.0;


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate heat transfer - update T
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar htc = 0.0;
    scalar T1 = calcHeatTransfer(td, dt, celli, htc, dhTrans);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate velocity - update U
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar Cud = 0.0;
    const vector U1 = calcVelocity(td, dt, Cud, dUTrans);


    // ~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate mass transfer
    // ~~~~~~~~~~~~~~~~~~~~~~~
    calcMassTransfer(td, dt, T0, T1, dMassMT);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate surface reactions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    calcSurfaceReactions(td, dt, celli, T0, T1, dMassMTSR, dMassSR);

    // New total mass
    const scalar mass1 = mass0 - sum(dMassMT) - dMassMTSR;

    // Ratio of mass devolatilised to the total volatile mass of the particle
    const scalar fVol = 1 -
        (YMixture_[0]*mass1)
       /(td.cloud().composition().YMixture0()[0]*mass0_);

    // Specific heat capacity of non-volatile components
    const scalar cpNonVolatile =
        (
            YMixture_[1]*td.cloud().composition().cpLiquid(YLiquid_, pc_, this->Tc_)
          + YMixture_[2]*td.cloud().composition().cpSolid(YSolid_)
        )/(YMixture_[1] + YMixture_[2]);

    // New specific heat capacity - linear variation until volatiles
    // have evolved
    const scalar cp1 = (cpNonVolatile - td.constProps().cp0())*fVol
       + td.constProps().cp0();


    // ~~~~~~~~~~~~~~~~~~~~~~~
    // Accumulate source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~

    // Transfer mass lost from particle to carrier mass source
    forAll(dMassMT, i)
    {
        td.cloud().rhoTrans(i)[celli] += np0*(dMassMT[i] + dMassSR[i]);
    }

    // Update momentum transfer
    td.cloud().UTrans()[celli] += np0*dUTrans;

    // Accumulate coefficient to be applied in carrier phase momentum coupling
    td.cloud().UCoeff()[celli] += np0*mass0*Cud;

    // Update enthalpy transfer
//    td.cloud().hTrans()[celli] += np0*(mass0*cp0*T0 - mass1*cp1*T1);
    td.cloud().hTrans()[celli] += np0*((mass0*cp0 - mass1*cp1)*T0 + dhTrans);

    // Accumulate coefficient to be applied in carrier phase enthalpy coupling
    td.cloud().hCoeff()[celli] += np0*htc*this->areaS();


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Remove the particle when mass falls below minimum threshold
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (mass1 < td.constProps().minParticleMass())
    {
        td.keepParticle = false;

        // Absorb particle(s) into carrier phase
        forAll(dMassMT, i)
        {
            td.cloud().rhoTrans(i)[celli] += np0*dMassMT[i];
        }
        td.cloud().hTrans()[celli] += np0*mass1*cp1*T1;
        td.cloud().UTrans()[celli] += np0*mass1*U1;
    }
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else
    {
        this->U_ = U1;
        this->T_ = T1;
        this->cp_ = cp1;

        // Update particle density or diameter
        if (td.cloud().massTransfer().changesVolume())
        {
            this->d_ = cbrt(mass1/this->rho_*6.0/mathematicalConstant::pi);
        }
        else
        {
            this->rho_ = mass1/this->volume();
        }
    }
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::calcUncoupled
(
    TrackData& td,
    const scalar dt,
    const label celli
)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Define local properties at beginning of timestep
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar T0 = this->T_;
    const scalar mass0 = this->mass();
//    const scalar cp0 = this->cp();

    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // Initialise transfer terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = vector::zero;

    // Enthalpy transfer from the particle to the carrier phase
    scalar dhTrans = 0.0;

    // Mass transfer from particle to carrier phase
    // - components exist in particle already
    scalarList dMassMT(td.cloud().gases().size(), 0.0);

    // Mass transfer due to surface reactions from particle to carrier phase
    // - components do not necessarily exist in particle already
    scalarList dMassSR(td.cloud().gases().size(), 0.0);

    // Total mass lost from particle due to surface reactions
    // - sub-model will adjust component mass fractions
    scalar dMassMTSR = 0.0;


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate heat transfer - update T
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar htc = 0.0;
    scalar T1 = calcHeatTransfer(td, dt, celli, htc, dhTrans);

    // Limit new temp max by vapourisarion temperature
    T1 = min(td.constProps().Tvap(), T1);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate velocity - update U
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar Cud = 0.0;
    const vector U1 = calcVelocity(td, dt, Cud, dUTrans);


    // ~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate mass transfer
    // ~~~~~~~~~~~~~~~~~~~~~~~
    calcMassTransfer(td, dt, T0, T1, dMassMT);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate surface reactions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    calcSurfaceReactions
    (
        td,
        dt,
        celli,
        T0,
        T1,
        dMassMTSR,
        dMassSR
    );

    // New total mass
    const scalar mass1 = mass0 - sum(dMassMT) - dMassMTSR;

    // Ratio of mass devolatilised to the total volatile mass of the particle
    const scalar fVol = 1 -
        (YMixture_[0]*mass1)
       /(td.cloud().composition().YMixture0()[0]*mass0_);

    // Specific heat capacity of non-volatile components
    const scalar cpNonVolatile =
        (
            YMixture_[1]*td.cloud().composition().cpLiquid(YLiquid_, pc_, this->Tc_)
          + YMixture_[2]*td.cloud().composition().cpSolid(YSolid_)
        )/(YMixture_[1] + YMixture_[2]);

    // New specific heat capacity - linear variation until volatiles
    // have evolved
    const scalar cp1 = (cpNonVolatile - td.constProps().cp0())*fVol
       + td.constProps().cp0();


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Remove the particle when mass falls below minimum threshold
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (mass1 < td.constProps().minParticleMass())
    {
        td.keepParticle = false;
    }
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else
    {
        this->U_ = U1;
        this->T_ = T1;
        this->cp_ = cp1;

        // Update particle density or diameter
        if (td.cloud().massTransfer().changesVolume())
        {
            this->d_ = cbrt(mass1/this->rho_*6.0/mathematicalConstant::pi);
        }
        else
        {
            this->rho_ = mass1/this->volume();
        }
    }
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::calcMassTransfer
(
    TrackData& td,
    const scalar dt,
    const scalar T0,
    const scalar T1,
    scalarList& dMassMT
)
{
    if (td.cloud().composition().YMixture0()[1]>SMALL)
    {
        notImplemented
        (
            "void Foam::ReactingParcel<ParcelType>::"
            "calcMassTransfer(...): no treatment currently "
            "available for particles containing liquid species"
        )
    }

    // Check that model is active, and that the parcel temperature is
    // within necessary limits for mass transfer to occur
    if
    (
        !td.cloud().massTransfer().active()
     || this->T_<td.constProps().Tvap()
     || this->T_<td.constProps().Tbp()
    )
    {
        return;
    }

    // Determine mass to add to carrier phase
    const scalar mass = this->mass();
    const scalar dMassTot = td.cloud().massTransfer().calculate
    (
        dt,
        mass0_,
        mass,
        td.cloud().composition().YMixture0(),
        YMixture_,
        T0,
        canCombust_
    );

    // Update (total) mass fractions
    YMixture_[0] = (YMixture_[0]*mass - dMassTot)/(mass - dMassTot);
    YMixture_[1] = YMixture_[1]*mass/(mass - dMassTot);
    YMixture_[2] = 1.0 - YMixture_[0] - YMixture_[1];

    // Add to cummulative mass transfer
    forAll (YGas_, i)
    {
        label id = td.cloud().composition().gasGlobalIds()[i];

        // Mass transfer
        scalar volatileMass = YGas_[i]*dMassTot;
        dMassMT[id] += volatileMass;
    }
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::calcSurfaceReactions
(
    TrackData& td,
    const scalar dt,
    const label celli,
    const scalar T0,
    const scalar T1,
    scalar& dMassMTSR,
    scalarList& dMassMT
)
{
    // Check that model is active
    if (!td.cloud().surfaceReaction().active() || !canCombust_)
    {
        return;
    }

    // Update mass transfer(s)
    // - Also updates Y()'s
    td.cloud().surfaceReaction().calculate
    (
        dt,
        celli,
        this->d_,
        T0,
        T1,
        this->Tc_,
        this->rhoc_,
        this->mass(),
        YGas_,
        YLiquid_,
        YSolid_,
        YMixture_,
        dMassMTSR,
        dMassMT
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "ReactingParcelIO.C"

// ************************************************************************* //

