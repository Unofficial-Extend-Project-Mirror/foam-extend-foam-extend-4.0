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

#include "ReactingCloud.H"
#include "CompositionModel.H"
#include "MassTransferModel.H"
#include "SurfaceReactionModel.H"

#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingCloud<ParcelType>::ReactingCloud
(
    const word& cloudType,
    const volPointInterpolation& vpi,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    hCombustionThermo& thermo,
    PtrList<specieReactingProperties>& gases
)
:
    ThermoCloud<ParcelType>(cloudType, vpi, rho, U, g, thermo),
    reactingCloud(),
    constProps_(this->particleProperties()),
    carrierThermo_(thermo),
    gases_(gases),
    compositionModel_
    (
        CompositionModel<ReactingCloud<ParcelType> >::New
        (
            this->particleProperties(),
            *this
        )
    ),
    massTransferModel_
    (
        MassTransferModel<ReactingCloud<ParcelType> >::New
        (
            this->particleProperties(),
            *this
        )
    ),
    surfaceReactionModel_
    (
        SurfaceReactionModel<ReactingCloud<ParcelType> >::New
        (
            this->particleProperties(),
            *this
        )
    ),
    rhoTrans_(thermo.composition().Y().size())
{
    // Set storage for mass source fields and initialise to zero
    forAll(rhoTrans_, i)
    {
        rhoTrans_.set
        (
            i,
            new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    this->name() + "rhoTrans" + Foam::name(i),
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh(),
                dimensionedScalar("zero", dimMass, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingCloud<ParcelType>::~ReactingCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::resetSourceTerms()
{
    ThermoCloud<ParcelType>::resetSourceTerms();
    forAll(rhoTrans_, i)
    {
        rhoTrans_[i].field() = 0.0;
    }
}


template<class ParcelType>
void Foam::ReactingCloud<ParcelType>::evolve()
{
    const volScalarField& T = carrierThermo_.T();
    const volScalarField cp = carrierThermo_.Cp();
    const volScalarField& p = carrierThermo_.p();

    autoPtr<interpolation<scalar> > rhoInterpolator = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        this->vpi(),
        this->rho()
    );

    autoPtr<interpolation<vector> > UInterpolator = interpolation<vector>::New
    (
        this->interpolationSchemes(),
        this->vpi(),
        this->U()
    );

    autoPtr<interpolation<scalar> > muInterpolator = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        this->vpi(),
        this->mu()
    );

    autoPtr<interpolation<scalar> > TInterpolator = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        this->vpi(),
        T
    );

    autoPtr<interpolation<scalar> > cpInterpolator = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        this->vpi(),
        cp
    );

    autoPtr<interpolation<scalar> > pInterpolator = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        this->vpi(),
        p
    );

    typename ParcelType::trackData td
    (
        *this,
        constProps_,
        rhoInterpolator(),
        UInterpolator(),
        muInterpolator(),
        TInterpolator(),
        cpInterpolator(),
        pInterpolator(),
        this->g().value()
    );

    inject(td);

    this->move(td);
}


template<class ParcelType>
template<class TrackingData>
void Foam::ReactingCloud<ParcelType>::inject
(
    TrackingData& td
)
{
    scalar time = this->db().time().value();

    scalar pRho = td.constProps().rho0();

    this->injection().prepareForNextTimeStep(this->time0(), time);

    // Number of parcels to introduce during this timestep
    const label nParcels = this->injection().nParcels();

    // Return if no parcels are required
    if (!nParcels)
    {
        this->postInjectCheck();
        return;
    }

    // Volume of particles to introduce during this timestep
    scalar pVolume = this->injection().volume();

    // Volume fraction to introduce during this timestep
    scalar pVolumeFraction = this->injection().volumeFraction();

    // Duration of injection period during this timestep
    scalar deltaT = min
    (
        this->db().time().deltaT().value(),
        min
        (
            time - this->injection().timeStart(),
            this->injection().timeEnd() - this->time0()
        )
    );

    // Pad injection time if injection starts during this timestep
    scalar padTime = max
    (
        0.0,
        this->injection().timeStart() - this->time0()
    );

    // Introduce new parcels linearly with time
    for (label iParcel=0; iParcel<nParcels; iParcel++)
    {
        // Calculate the pseudo time of injection for parcel 'iParcel'
        scalar timeInj = this->time0() + padTime + deltaT*iParcel/nParcels;

        // Determine injected parcel properties
        vector pPosition = this->injection().position
        (
            iParcel,
            timeInj,
            this->meshInfo()
        );

        // Diameter of parcels
        scalar pDiameter = this->injection().d0(iParcel, timeInj);

        // Number of particles per parcel
        scalar pNumberOfParticles = this->setNumberOfParticles
        (
            nParcels,
            pDiameter,
            pVolumeFraction,
            pRho,
            pVolume
        );

        // Velocity of parcels
        vector pU = this->injection().velocity
        (
            iParcel,
            timeInj,
            this->meshInfo()
        );

        // Determine the injection cell
        label pCell = -1;
        this->setInjectorCellAndPosition(pCell, pPosition);

        if (pCell >= 0)
        {
            // construct the parcel that is to be injected
            ParcelType* pPtr = new ParcelType
            (
                td.cloud(),
                this->parcelTypeId(),
                pPosition,
                pCell,
                pDiameter,
                pU,
                pNumberOfParticles,
                composition().YGas0(),
                composition().YLiquid0(),
                composition().YSolid0(),
                composition().YMixture0(),
                td.constProps()
            );

            scalar dt = time - timeInj;

            pPtr->stepFraction() = (this->db().time().deltaT().value() - dt)
                /this->db().time().deltaT().value();

            this->injectParcel(td, pPtr);
         }
    }

    this->postInjectCheck();

    if (debug)
    {
        this->dumpParticlePositions();
    }
}


// ************************************************************************* //
