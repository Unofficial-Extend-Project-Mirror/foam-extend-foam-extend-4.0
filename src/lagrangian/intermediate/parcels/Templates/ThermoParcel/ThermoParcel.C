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

#include "ThermoParcel.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::ThermoParcel<ParcelType>::updateCellQuantities
(
    TrackData& td,
    const scalar dt,
    const label celli
)
{
    KinematicParcel<ParcelType>::updateCellQuantities(td, dt, celli);

    Tc_ = td.TInterp().interpolate(this->position(), celli);
    cpc_ = td.cpInterp().interpolate(this->position(), celli);
}


template<class ParcelType>
template<class TrackData>
void Foam::ThermoParcel<ParcelType>::calcCoupled
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
//    const scalar T0 = T_;
//    const scalar cp0 = cp_;


    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // Initialise transfer terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = vector::zero;

    // Enthalpy transfer from the particle to the carrier phase
    scalar dhTrans = 0.0;


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate velocity - update U
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar Cud = 0.0;
    const vector U1 = calcVelocity(td, dt, Cud, dUTrans);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate heat transfer - update T
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar htc = 0.0;
    const scalar T1 = calcHeatTransfer(td, dt, celli, htc, dhTrans);


    // ~~~~~~~~~~~~~~~~~~~~~~~
    // Accumulate source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~

    // Update momentum transfer
    td.cloud().UTrans()[celli] += np0*dUTrans;

    // Accumulate coefficient to be applied in carrier phase momentum coupling
    td.cloud().UCoeff()[celli] += np0*mass0*Cud;

    // Update enthalpy transfer
    td.cloud().hTrans()[celli] += np0*dhTrans;

    // Accumulate coefficient to be applied in carrier phase enthalpy coupling
    td.cloud().hCoeff()[celli] += np0*htc*this->areaS();


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    this->U() = U1;
    this->T() = T1;
}


template<class ParcelType>
template<class TrackData>
void Foam::ThermoParcel<ParcelType>::calcUncoupled
(
    TrackData& td,
    const scalar dt,
    const label celli
)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // Initialise transfer terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = vector::zero;

    // Enthalpy transfer from the particle to the carrier phase
    scalar dhTrans = 0.0;


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate velocity - update U
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar Cud = 0.0;
    this->U_ = calcVelocity(td, dt, Cud, dUTrans);


    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate heat transfer - update T
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalar htc = 0.0;
    T_ = calcHeatTransfer(td, dt, celli, htc, dhTrans);
}


template<class ParcelType>
template <class TrackData>
Foam::scalar Foam::ThermoParcel<ParcelType>::calcHeatTransfer
(
    TrackData& td,
    const scalar dt,
    const label celli,
    scalar& htc,
    scalar& dhTrans
)
{
    if (!td.cloud().heatTransfer().active())
    {
        htc = 0.0;
        dhTrans = 0.0;
        return T_;
    }

    // Calc heat transfer coefficient
    htc = td.cloud().heatTransfer().h
    (
        this->d_,
        this->U_ - this->Uc_,
        this->rhoc_,
        this->rho_,
        cpc_,
        cp_,
        this->muc_
    );

    // Determine ap and bp coefficients
    scalar ap = Tc_;
    scalar bp = htc;
    if (td.cloud().radiation())
    {
        // Carrier phase incident radiation field
        // - The G field is not interpolated to the parcel position
        //   Instead, the cell centre value is applied directly
        const scalarField& G = td.cloud().mesh().objectRegistry
            ::lookupObject<volScalarField>("G");

        // Helper variables
        const scalar sigma = radiation::sigmaSB.value();
        const scalar epsilon = td.constProps().epsilon0();
        const scalar epsilonSigmaT3 = epsilon*sigma*pow3(T_);
        ap = (htc*Tc_ + 0.25*epsilon*G[celli])/(htc + epsilonSigmaT3);
        bp += epsilonSigmaT3;
    }
    bp *= 6.0/(this->rho_*this->d_*cp_);


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Set new particle temperature
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalar Tnew = td.cloud().TIntegrator().integrate(T_, dt, ap, bp);

    dhTrans = -this->mass()*cp_*(Tnew - T_);

    return Tnew;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "ThermoParcelIO.C"

// ************************************************************************* //

