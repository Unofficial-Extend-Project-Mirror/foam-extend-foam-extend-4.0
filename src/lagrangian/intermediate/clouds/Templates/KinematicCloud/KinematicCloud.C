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

#include "KinematicCloud.H"
#include "DispersionModel.H"
#include "DragModel.H"
#include "InjectionModel.H"
#include "WallInteractionModel.H"

#include "IntegrationScheme.H"

#include "interpolationCellPoint.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::setInjectorCellAndPosition
(
    label& pCell,
    vector& pPosition
)
{
    const vector originalPosition = pPosition;

    bool foundCell = false;

    pCell = mesh_.findCell(pPosition);

    if (pCell >= 0)
    {
        const vector& C = mesh_.C()[pCell];
        pPosition += 1.0e-6*(C - pPosition);

        foundCell = mesh_.pointInCell
        (
            pPosition,
            pCell
        );
    }
    reduce(foundCell, orOp<bool>());

    // Last chance - find nearest cell and try that one
    // - the point is probably on an edge
    if (!foundCell)
    {
        pCell =  mesh_.findNearestCell(pPosition);

        if (pCell >= 0)
        {
            const vector& C = mesh_.C()[pCell];
            pPosition += 1.0e-6*(C - pPosition);

            foundCell = mesh_.pointInCell
            (
                pPosition,
                pCell
            );
        }
        reduce(foundCell, orOp<bool>());
    }

    if (!foundCell)
    {
        FatalErrorIn
        (
            "void KinematicCloud<ParcelType>::findInjectorCell"
            "(label&, vector&)"
        )<< "Cannot find parcel injection cell. "
         << "Parcel position = " << originalPosition << nl
         << abort(FatalError);
    }
}


template<class ParcelType>
Foam::scalar Foam::KinematicCloud<ParcelType>::setNumberOfParticles
(
    const label nParcels,
    const scalar pDiameter,
    const scalar pVolumeFraction,
    const scalar pRho,
    const scalar pVolume
)
{
    scalar nP = 0.0;
    switch (parcelBasis_)
    {
        case pbMass:
        {
            nP = pVolumeFraction*massTotal_/nParcels
               /(pRho*mathematicalConstant::pi/6.0*pow(pDiameter, 3));
            break;
        }
        case pbNumber:
        {
            nP = pVolumeFraction*massTotal_/(pRho*pVolume);
            break;
        }
        default:
        {
            nP = 0.0;
            FatalErrorIn
            (
                "Foam::KinematicCloud<ParcelType>::setNumberOfParticles"
                "(const label, const scalar, const scalar, const scalar, "
                "const scalar)"
            )<< "Unknown parcelBasis type" << nl
             << exit(FatalError);
        }
    }

    return nP;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicCloud<ParcelType>::KinematicCloud
(
    const word& cloudType,
    const volPointInterpolation& vpi,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g
)
:
    Cloud<ParcelType>(rho.mesh(), cloudType, false),
    kinematicCloud(),
    cloudType_(cloudType),
    mesh_(rho.mesh()),
    vpi_(vpi),
    particleProperties_
    (
        IOobject
        (
            cloudType + "Properties",
            rho.mesh().time().constant(),
            rho.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    constProps_(particleProperties_),
    parcelTypeId_(readLabel(particleProperties_.lookup("parcelTypeId"))),
    coupled_(particleProperties_.lookup("coupled")),
    rndGen_(label(0)),
    time0_(this->db().time().value()),
    parcelBasisType_(particleProperties_.lookup("parcelBasisType")),
    parcelBasis_(pbNumber),
    massTotal_
    (
        dimensionedScalar(particleProperties_.lookup("massTotal")).value()
    ),
    massInjected_(0.0),
    rho_(rho),
    U_(U),
    mu_(mu),
    g_(g),
    interpolationSchemes_(particleProperties_.subDict("interpolationSchemes")),
    dispersionModel_
    (
        DispersionModel<KinematicCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    ),
    dragModel_
    (
        DragModel<KinematicCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    ),
    injectionModel_
    (
        InjectionModel<KinematicCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    ),
    wallInteractionModel_
    (
        WallInteractionModel<KinematicCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    ),
    UIntegrator_
    (
        vectorIntegrationScheme::New
        (
            "U",
            particleProperties_.subDict("integrationSchemes")
        )
    ),
    nInjections_(0),
    nParcelsAdded_(0),
    nParcelsAddedTotal_(0),
    UTrans_
    (
        IOobject
        (
            this->name() + "UTrans",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimensionSet(1, 1, -1, 0, 0), vector::zero)
    ),
    UCoeff_
    (
        IOobject
        (
            this->name() + "UCoeff",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, 0, -1, 0, 0), 0.0)
    )
{
    if (parcelBasisType_ == "mass")
    {
        parcelBasis_ = pbMass;
    }
    else if (parcelBasisType_ == "number")
    {
        parcelBasis_ = pbNumber;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::KinematicCloud<ParcelType>::KinematicCloud"
            "(const word&, const volPointInterpolation&, const volScalarField&"
            ", const volVectorField&, const volScalarField&, const "
            "dimensionedVector&)"
        )<< "parcelBasisType must be either 'number' or 'mass'" << nl
         << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicCloud<ParcelType>::~KinematicCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::resetSourceTerms()
{
    UTrans_.field() = vector::zero;
    UCoeff_.field() = 0.0;
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::evolve()
{
    autoPtr<interpolation<scalar> > rhoInterpolator =
        interpolation<scalar>::New
        (
            interpolationSchemes_,
            vpi_,
            rho_
        );

    autoPtr<interpolation<vector> > UInterpolator =
        interpolation<vector>::New
        (
            interpolationSchemes_,
            vpi_,
            U_
        );

    autoPtr<interpolation<scalar> > muInterpolator =
        interpolation<scalar>::New
        (
            interpolationSchemes_,
            vpi_,
            mu_
        );

    typename ParcelType::trackData td
    (
        *this,
        constProps_,
        rhoInterpolator(),
        UInterpolator(),
        muInterpolator(),
        g_.value()
    );

    inject(td);

    move(td);
}


template<class ParcelType>
template<class TrackingData>
void Foam::KinematicCloud<ParcelType>::move
(
    TrackingData& td
)
{
    if (coupled_)
    {
        resetSourceTerms();
    }
    Cloud<ParcelType>::move(td);
}


template<class ParcelType>
template<class TrackingData>
void Foam::KinematicCloud<ParcelType>::inject
(
    TrackingData& td
)
{
    scalar time = this->db().time().value();

    scalar pRho = td.constProps().rho0();

    this->injection().prepareForNextTimeStep(time0_, time);

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
            this->injection().timeEnd() - time0_
        )
    );

    // Pad injection time if injection starts during this timestep
    scalar padTime = max
    (
        0.0,
        this->injection().timeStart() - time0_
    );

    // Introduce new parcels linearly with time
    for (label iParcel=0; iParcel<nParcels; iParcel++)
    {
        // Calculate the pseudo time of injection for parcel 'iParcel'
        scalar timeInj = time0_ + padTime + deltaT*iParcel/nParcels;

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
        scalar pNumberOfParticles = setNumberOfParticles
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
        setInjectorCellAndPosition(pCell, pPosition);

        if (pCell >= 0)
        {
            // construct the parcel that is to be injected
            ParcelType* pPtr = new ParcelType
            (
                td.cloud(),
                parcelTypeId_,
                pPosition,
                pCell,
                pDiameter,
                pU,
                pNumberOfParticles,
                td.constProps()
            );

            scalar dt = time - timeInj;

            pPtr->stepFraction() = (this->db().time().deltaT().value() - dt)
                /this->time().deltaT().value();

            this->injectParcel(td, pPtr);
         }
    }

    this->postInjectCheck();

    if (debug)
    {
        this->dumpParticlePositions();
    }
}


template<class ParcelType>
template<class TrackingData>
void Foam::KinematicCloud<ParcelType>::injectParcel
(
    TrackingData& td,
    ParcelType* p
)
{
    addParticle(p);
    nParcelsAdded_++;
    nParcelsAddedTotal_++;
    massInjected_ += p->mass()*p->nParticle();
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::postInjectCheck()
{
    if (nParcelsAdded_)
    {
        Pout<< "\n--> Cloud: " << this->name() << nl
            << "    Added " << nParcelsAdded_
            <<  " new parcels" << nl << endl;
    }

    // Reset parcel counters
    nParcelsAdded_ = 0;

    // Set time for start of next injection
    time0_ = this->db().time().value();

    // Increment number of injections
    nInjections_++;
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::info() const
{
    Info<< "Cloud name: " << this->name() << nl
        << "    Parcels added during this run   = "
        << returnReduce(nParcelsAddedTotal_, sumOp<label>()) << nl
        << "    Mass introduced during this run = "
        << returnReduce(massInjected_, sumOp<scalar>()) << nl
        << "    Current number of parcels       = "
        << returnReduce(this->size(), sumOp<label>()) << nl
        << "    Current mass in system          = "
        << returnReduce(massInSystem(), sumOp<scalar>()) << nl
        << endl;
}


template<class ParcelType>
void Foam::KinematicCloud<ParcelType>::dumpParticlePositions() const
{
    OFstream pObj
    (
        this->db().time().path()/"parcelPositions_"
      + this->name() + "_"
      + name(this->nInjections_) + ".obj"
    );

    forAllConstIter(typename KinematicCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();
        pObj<< "v " << p.position().x() << " " << p.position().y() << " "
            << p.position().z() << nl;
    }

    pObj.flush();
}


// ************************************************************************* //
