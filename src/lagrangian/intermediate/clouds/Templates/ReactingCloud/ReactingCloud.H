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
    Foam::ReactingCloud

Description
    Templated base class for reactive cloud
    - Adds to kinematic cloud
      - Heat transfer

SourceFiles
    ReactingCloudI.H
    ReactingCloud.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingCloud_H
#define ReactingCloud_H

#include "autoPtr.H"
#include "hCombustionThermo.H"

#include "ThermoCloud.H"
#include "reactingCloud.H"

#include "ReactingCloudThermoTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class CompositionModel;

template<class CloudType>
class MassTransferModel;

template<class CloudType>
class SurfaceReactionModel;

/*---------------------------------------------------------------------------*\
                        Class ReactingCloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ReactingCloud
:
    public ThermoCloud<ParcelType>,
    public reactingCloud
{

    // Private data

        //- Parcel constant properties
        typename ParcelType::constantProperties constProps_;

        //- Thermodynamics package (combustion)
        hCombustionThermo& carrierThermo_;

        //- Gas phase properties
        PtrList<specieReactingProperties>& gases_;

        // References to the cloud sub-models

            //- Reacting composition model
            autoPtr<CompositionModel<ReactingCloud<ParcelType> > >
                compositionModel_;

            //- Mass transfer model
            autoPtr<MassTransferModel<ReactingCloud<ParcelType> > >
                massTransferModel_;

            //- Surface reaction model
            autoPtr<SurfaceReactionModel<ReactingCloud<ParcelType> > >
                surfaceReactionModel_;

        // Sources

            //- Mass transfer fields - one per carrier phase specie
            PtrList<DimensionedField<scalar, volMesh> > rhoTrans_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ReactingCloud(const ReactingCloud&);

        //- Disallow default bitwise assignment
        void operator=(const ReactingCloud&);


protected:

    //- Inject more parcels
    template<class TrackingData>
    void inject(TrackingData& td);


public:

    //- Runtime type information
//    TypeName("ReactingCloud");


    // Constructors

        //- Construct given carrier gas fields
        ReactingCloud
        (
            const word& cloudType,
            const volPointInterpolation& vpi,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            hCombustionThermo& thermo,
            PtrList<specieReactingProperties>& gases
        );


    // Destructor

        virtual ~ReactingCloud();


    // Member Functions

        // Access

            //- Return const access to carrier phase thermo package
            inline const hCombustionThermo& carrierThermo() const;

            //- Return access to carrier phase thermo package
            inline hCombustionThermo& carrierThermo();

            //- Gas phase properties
            inline const PtrList<specieReactingProperties>& gases() const;


            // Sub-models

                //- Return reference to reacting composition model
                inline const CompositionModel<ReactingCloud<ParcelType> >&
                    composition() const;

                //- Return reference to mass transfer model
                inline const MassTransferModel<ReactingCloud<ParcelType> >&
                    massTransfer() const;

                //- Return reference to surface reaction model
                inline const SurfaceReactionModel<ReactingCloud<ParcelType> >&
                    surfaceReaction() const;


            // Sources

                //- Mass

                    //- Return reference to mass source for field i
                    inline DimensionedField<scalar, volMesh>&
                        rhoTrans(const label i);

                    //- Return reference to mass source fields
                    inline PtrList<DimensionedField<scalar, volMesh> >&
                        rhoTrans();

                    //- Return tmp mass source for field i
                    //  Fully explicit
                    inline tmp<DimensionedField<scalar, volMesh> >
                        Srho1(const label i) const;

                    //- Return tmp total mass source for carrier phase
                    //  Fully explicit
                    inline tmp<DimensionedField<scalar, volMesh> >
                        Srho1() const;


        // Cloud evolution functions

            //- Reset the spray source terms
            void resetSourceTerms();

            //- Evolve the spray (move, inject)
            void evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ReactingCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
