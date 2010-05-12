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

Class
    Foam::KinematicCloud

Description
    Templated base class for kinematic cloud

    - Kinematic only
      - Dispersion model
      - Drag model
      - Injection model
      - Wall interaction model

SourceFiles
    KinematicCloudI.H
    KinematicCloud.C

\*---------------------------------------------------------------------------*/

#ifndef KinematicCloud_H
#define KinematicCloud_H

#include "Cloud.H"
#include "kinematicCloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "Random.H"
#include "fvMesh.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "fvMatrices.H"
#include "fvm.H"

#include "IntegrationSchemesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class DispersionModel;

template<class CloudType>
class DragModel;

template<class CloudType>
class InjectionModel;

template<class CloudType>
class WallInteractionModel;

/*---------------------------------------------------------------------------*\
                       Class KinematicCloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class KinematicCloud
:
    public Cloud<ParcelType>,
    public kinematicCloud
{

public:

    // Enumerations

        //- Parcel basis representation options
        //  i.e constant number of particles OR constant mass per parcel
        enum parcelBasis
        {
            pbNumber,
            pbMass
        };


private:

    // Private data

        //- Cloud type - used to set the name of the parcel properties
        //  dictionary by appending "Properties"
        const word cloudType_;

        //- References to the mesh and time databases
        const fvMesh& mesh_;
//        const Time& runTime_;

        //- Reference to the interpolation for the carrier phase to the parcels
        const volPointInterpolation& vpi_;

        //- Dictionary of particle properties
        IOdictionary particleProperties_;

        //- Parcel constant properties
        typename ParcelType::constantProperties constProps_;

        //- Parcel type id - used to flag the type of parcels issued by this
        //  cloud
        const label parcelTypeId_;

        //- Flag to indicate whether parcels are coupled to the carrier phase
        //  i.e. whether or not to generate source terms for carrier phase
        Switch coupled_;

        //- Random number generator - used by some injection routines
        Random rndGen_;

        //- Time at beginning of timestep
        scalar time0_;


        // Injection properties

            //- Parcel basis
            const word parcelBasisType_;
            parcelBasis parcelBasis_;

            //- Total mass to inject [kg]
            scalar massTotal_;

            //- Total mass injected to date [kg]
            scalar massInjected_;


        // References to the carrier gas fields

            //- Density
            const volScalarField& rho_;

            //- Velocity
            const volVectorField& U_;

            //- Dynamic viscosity
            const volScalarField& mu_;


        // Environmental properties

            //- Gravity
            const dimensionedVector& g_;


        // Interpolation

            dictionary interpolationSchemes_;


        // References to the cloud sub-models

            //- Dispersion model
            autoPtr<DispersionModel<KinematicCloud<ParcelType> > >
                dispersionModel_;

            //- Drag transfer model
            autoPtr<DragModel<KinematicCloud<ParcelType> > > dragModel_;

            //- Injector model
            autoPtr<InjectionModel<KinematicCloud<ParcelType> > >
                injectionModel_;

            //- Wall interaction model
            autoPtr<WallInteractionModel<KinematicCloud<ParcelType> > >
                wallInteractionModel_;


        // Reference to the particle integration schemes

            //- Velocity integration
            autoPtr<vectorIntegrationScheme> UIntegrator_;


        // Counters

            //- Number of injections counter
            label nInjections_;

            //- Running counters describing parcels added during each
            //  injection
            label nParcelsAdded_;
            label nParcelsAddedTotal_;


        // Sources

            //- Momentum
            DimensionedField<vector, volMesh> UTrans_;

            //- Coefficient for carrier phase U equation
            DimensionedField<scalar, volMesh> UCoeff_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        KinematicCloud(const KinematicCloud&);

        //- Disallow default bitwise assignment
        void operator=(const KinematicCloud&);


protected:

    // Protected member functions

        //- Set parcel position and cell into which parcel is introduced
        void setInjectorCellAndPosition
        (
            label& pCell,
            vector& pPosition
        );

        //- Set the number of particles per parcel
        scalar setNumberOfParticles
        (
            const label nParcels,
            const scalar pDiameter,
            const scalar pVolumeFraction,
            const scalar pRho,
            const scalar pVolume
        );

        //- Inject more parcels
        template<class TrackingData>
        void inject(TrackingData& td);

        //- Inject parcel if it is valid - delete otherwise
        template<class TrackingData>
        void injectParcel
        (
            TrackingData& td,
            ParcelType* p
        );

        //- Move the parcels
        template<class TrackingData>
        void move(TrackingData& td);

        //- Post-injection checks
        void postInjectCheck();


public:

    //- Runtime type information
//    TypeName("KinematicCloud");

    // Constructors

        //- Construct given carrier gas fields
        KinematicCloud
        (
            const word& cloudType,
            const volPointInterpolation& vpi,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g
        );


    // Destructor

        virtual ~KinematicCloud();


    // Member Functions

        // Access

            // References to the mesh and databases

                //- Return the cloud type
                inline const word& cloudType() const;

                //- Return the parcel type id
                inline label parcelTypeId() const;

                //- Return refernce to the mesh
                inline const fvMesh& mesh() const;

                //- Retuen reference to the interpolation
                inline const volPointInterpolation& vpi() const;

                //- Return particle properties dictionary
                inline const IOdictionary& particleProperties() const;

            //- Return coupled flag
            inline const Switch coupled() const;

            //- Return refernce to the random object
            inline Random& rndGen();

            //- Return the start of injection interval time
            inline scalar time0() const;

            //- Return a reference to the mass of particles to introduce
            inline scalar massTotal() const;


            // References to the carrier gas fields

                //- Return carrier gas velocity
                inline const volVectorField& U() const;

                //- Return carrier gas density
                inline const volScalarField& rho() const;

                //- Return carrier gas dynamic viscosity
                inline const volScalarField& mu() const;


            // Environmental properties

                //- Gravity
                inline const dimensionedVector& g() const;


            // Interpolations

                //- Return reference to the interpolation dictionary
                inline const dictionary& interpolationSchemes() const;


            // Sub-models

                //- Return reference to dispersion model
                inline const DispersionModel<KinematicCloud<ParcelType> >&
                    dispersion() const;

                inline DispersionModel<KinematicCloud<ParcelType> >&
                    dispersion();

                //- Return reference to drag model
                inline const DragModel<KinematicCloud<ParcelType> >&
                    drag() const;

                //- Return reference to injection model
                inline const InjectionModel<KinematicCloud<ParcelType> >&
                    injection() const;

                inline InjectionModel<KinematicCloud<ParcelType> >&
                    injection();

                //- Return reference to wall interaction model
                inline const WallInteractionModel<KinematicCloud<ParcelType> >&
                    wallInteraction() const;


            // Integration schemes

                //-Return reference to velocity integration
                inline const vectorIntegrationScheme& UIntegrator() const;


            // Sources

                // Momentum

                    //- Return reference to momentum source
                    inline DimensionedField<vector, volMesh>& UTrans();

                    //- Coefficient for carrier phase U equation
                    inline DimensionedField<scalar, volMesh>& UCoeff();

                    //- Return tmp momentum source term - fully explicit
                    inline tmp<DimensionedField<vector, volMesh> > SU1() const;

                    //- Return tmp momentum source term - semi-implicit
                    inline tmp<fvVectorMatrix> SU2(volVectorField& U) const;



        // Check

            //- Total mass injected
            inline scalar massInjected() const;

            //- Total mass in system
            inline scalar massInSystem() const;

            //- Print cloud information
            void info() const;

            //- Dump particle positions to .obj file
            void dumpParticlePositions() const;


            // Counters

                //- Return the number of injections
                inline label nInjections() const;

                //- Return the total number parcels added
                inline label nParcelsAddedTotal() const;


            // Fields

                //- Return the particle volume fraction field
                //  Note: for particles belonging to this cloud only
                inline const tmp<volScalarField> theta() const;

                //- Return the particle mass fraction field
                //  Note: for particles belonging to this cloud only
                inline const tmp<volScalarField> alpha() const;


        // Cloud evolution functions

            //- Reset the spray source terms
            void resetSourceTerms();

            //- Evolve the spray (move, inject)
            void evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "KinematicCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "KinematicCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
