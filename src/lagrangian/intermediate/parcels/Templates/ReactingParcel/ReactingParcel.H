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
    Foam::ReactingParcel

Description
    Reactinf parcel class with one/two-way coupling with the continuous
    phase. Includes thermo parcel sub-models, plus:
    - combustion

SourceFiles
    ReactingParcelI.H
    ReactingParcel.C
    ReactingParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingParcel_H
#define ReactingParcel_H

#include "IOstream.H"
#include "autoPtr.H"
#include "interpolationCellPoint.H"
#include "contiguous.H"

#include "ThermoParcel.H"
#include "ReactingCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class ReactingParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const ReactingParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                        Class ReactingParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ReactingParcel
:
    public ThermoParcel<ParcelType>
{

public:

    //- Class to hold reacting particle constant properties
    class constantProperties
    :
        public ThermoParcel<ParcelType>::constantProperties
    {

        // Private data

            //- Vapourisation temperature [K]
            const scalar Tvap_;

            //- Boiling point [K]
            const scalar Tbp_;


    public:

        //- Constructor
        constantProperties(const dictionary& dict);

        // Access

            //- Return const access to the vapourisation temperature
            inline scalar Tvap() const;

            //- Return const access to the boiling point
            inline scalar Tbp() const;
    };


    //- Class used to pass reacting tracking data to the trackToFace function
    class trackData
    :
        public ThermoParcel<ParcelType>::trackData
    {

        // Private data

            //- Reference to the cloud containing this particle
            ReactingCloud<ParcelType>& cloud_;

            //- Particle constant properties
            const constantProperties& constProps_;

            //- Interpolator for continuous phase pressure field
            const interpolation<scalar>& pInterp_;


    public:

        // Constructors

            //- Construct from components
            inline trackData
            (
                ReactingCloud<ParcelType>& cloud,
                const constantProperties& constProps,
                const interpolation<scalar>& rhoInterp,
                const interpolation<vector>& UInterp,
                const interpolation<scalar>& muInterp,
                const interpolation<scalar>& TInterp,
                const interpolation<scalar>& CpInterp,
                const interpolation<scalar>& pInterp,
                const vector& g
            );


        // Member functions

            //- Return access to the owner cloud
            inline ReactingCloud<ParcelType>& cloud();

            //- Return const access to the constant properties
            inline const constantProperties& constProps() const;

            //- Return conat access to the interpolator for continuous
            //  phase pressure field
            inline const interpolation<scalar>& pInterp() const;
    };


protected:

    // Protected data

        // Parcel properties

            //- Initial particle mass [kg]
            scalar mass0_;

            //- Mass fractions of mixture []
            scalarField YMixture_;

            //- Mass fractions of gases []
            scalarField YGas_;

            //- Mass fractions of liquids []
            scalarField YLiquid_;

            //- Mass fractions of solids []
            scalarField YSolid_;

            //- Flag to say that the particle is allowed to combust
            //  Only true after volatile content falls below threshold value
            bool canCombust_;


        // Cell-based quantities

            //- Pressure [Pa]
            scalar pc_;


    // Protected member functions

        //- Calculate mass transfer
        template<class TrackData>
        void calcMassTransfer
        (
            TrackData& td,
            const scalar dt,
            const scalar T0,
            const scalar T1,
            scalarList& dMassMT
        );

        //- Calculate surface reactions
        template<class TrackData>
        void calcSurfaceReactions
        (
            TrackData& td,
            const scalar dt,
            const label celli,
            const scalar T0,
            const scalar T1,
            scalar& dMassMTSR,
            scalarList& dMassMT
        );


public:

    //- Runtime type information
    TypeName("ReactingParcel");

    friend class Cloud<ParcelType>;


    // Constructors

        //- Construct from components
        inline ReactingParcel
        (
            ReactingCloud<ParcelType>& owner,
            const label typeId,
            const vector& position,
            const label celli,
            const scalar d0,
            const vector& U0,
            const scalar nParticle0,
            const scalarField& YGas0,
            const scalarField& YLiquid0,
            const scalarField& YSolid0,
            const scalarField& YMixture0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        ReactingParcel
        (
            const Cloud<ParcelType>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<ParcelType> clone() const
        {
            return autoPtr<ParcelType>(new ReactingParcel<ParcelType>(*this));
        }


    // Member Functions

        // Access

            //- Return const access to mass fractions of gases
            inline const scalarField& YGas() const;

            //- Return const access to mass fractions of liquids
            inline const scalarField& YLiquid() const;

            //- Return const access to mass fractions of solids
            inline const scalarField& YSolid() const;

            //- Return const access to mass fractions of mixture
            inline const scalarField& YMixture() const;

            //- Return const access to initial mass
            inline scalar mass0() const;


        // Edit

            //- Return access to mass fractions of gases
            inline scalarField& YGas();

            //- Return access to mass fractions of liquids
            inline scalarField& YLiquid();

            //- Return access to mass fractions of solids
            inline scalarField& YSolid();

            //- Return access to mass fractions of mixture
            inline scalarField& YMixture();

            //- Return access to initial mass
            inline scalar& mass0();


        // Main calculation loop

            //- Update cell based quantities
            template<class TrackData>
            void updateCellQuantities
            (
                TrackData& td,
                const scalar dt,
                const label celli
            );

            //- Coupled calculation with the continuous phase
            template<class TrackData>
            void calcCoupled
            (
                TrackData& td,
                const scalar dt,
                const label celli
            );

            //- Uncoupled calculation with the continuous phase
            template<class TrackData>
            void calcUncoupled
            (
                TrackData& td,
                const scalar dt,
                const label
            );


        // I-O

            static void readFields(ReactingCloud<ParcelType>& c);

            static void writeFields(const ReactingCloud<ParcelType>& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const ReactingParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

