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
    Foam::ThermoCloud

Description
    Templated base class for thermodynamic cloud
    - Adds to kinematic cloud
      - Heat transfer

SourceFiles
    ThermoCloudI.H
    ThermoCloud.C

\*---------------------------------------------------------------------------*/

#ifndef ThermoCloud_H
#define ThermoCloud_H

#include "autoPtr.H"
#include "hCombustionThermo.H"

#include "KinematicCloud.H"
#include "thermoCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class HeatTransferModel;

/*---------------------------------------------------------------------------*\
                        Class ThermoCloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ThermoCloud
:
    public KinematicCloud<ParcelType>,
    public thermoCloud
{

    // Private data

        //- Thermo parcel constant properties
        typename ParcelType::constantProperties constProps_;


        // References to the carrier gas fields

            //- Thermodynamics package (basic)
            basicThermo& carrierThermo_;


        // References to the cloud sub-models

            //- Heat transfer model
            autoPtr<HeatTransferModel<ThermoCloud<ParcelType> > >
                heatTransferModel_;


        // Reference to the particle integration schemes

            //- Temperature integration
            autoPtr<scalarIntegrationScheme> TIntegrator_;


        // Modelling options

            //- Include radiation
            Switch radiation_;


        // Sources

            //- Enthalpy transfer
            DimensionedField<scalar, volMesh> hTrans_;

            //- Coefficient for carrier phase h equation
            DimensionedField<scalar, volMesh> hCoeff_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ThermoCloud(const ThermoCloud&);

        //- Disallow default bitwise assignment
        void operator=(const ThermoCloud&);


protected:

    //- Inject more parcels
    template<class TrackingData>
    void inject(TrackingData& td);


public:

    //- Runtime type information
//    TypeName("ThermoCloud");


    // Constructors

        //- Construct given carrier gas fields
        ThermoCloud
        (
            const word& cloudType,
            const volPointInterpolation& vpi,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            basicThermo& thermo
        );


    // Destructor

        virtual ~ThermoCloud();


    // Member Functions

        // Access

            //- Return const access to thermo package
            inline const basicThermo& carrierThermo() const;

            //- Return access to thermo package
            inline basicThermo& carrierThermo();


            // Sub-models

                //- Return reference to heat transfer model
                inline const HeatTransferModel<ThermoCloud<ParcelType> >&
                    heatTransfer() const;


            // Integration schemes

                //-Return reference to velocity integration
                inline const scalarIntegrationScheme& TIntegrator() const;


            // Modelling options

                //- Radiation flag
                inline bool radiation() const;


            // Sources

                // Enthalpy

                    //- Return reference to enthalpy source
                    inline DimensionedField<scalar, volMesh>& hTrans();

                    //- Coefficient for carrier phase h equation
                    inline DimensionedField<scalar, volMesh>& hCoeff();

                    //- return tmp enthalpy source term - fully explicit
                    inline tmp<DimensionedField<scalar, volMesh> > Sh1() const;

                    //- Return tmp enthalpy source term - semi-implicit
                    inline tmp<fvScalarMatrix> Sh2(volScalarField& h) const;


                // Radiation - overrides thermoCloud virtual abstract members

                    //- Return tmp equivalent particulate emission
                    inline tmp<volScalarField> Ep() const;

                    //- Return tmp equivalent particulate absorption
                    inline tmp<volScalarField> ap() const;

                    //- Return tmp equivalent particulate scattering factor
                    inline tmp<volScalarField> sigmap() const;


        // Cloud evolution functions

            //- Evolve the spray (move, inject)
            void evolve();

            //- Reset the spray source terms
            void resetSourceTerms();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ThermoCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ThermoCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
