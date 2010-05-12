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
    Foam::InjectionModel


Description
    Templated injection model class

SourceFiles
    InjectionModel.C
    NewInjectionModel.C

\*---------------------------------------------------------------------------*/

#ifndef InjectionModel_H
#define InjectionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class InjectionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class InjectionModel
{

    // Private data

        //- The cloud dictionary
        const dictionary& dict_;

        // Reference to the owner cloud class
        CloudType& owner_;

        //- The coefficients dictionary
        const dictionary coeffDict_;


protected:

    // Protected data

        // Global injection properties

            //- Start of injection [s]
            scalar SOI_;

            //- Total volume of parcels to introduce [m^3]
            //  Initialised in the individual injection models
            scalar volumeTotal_;


        // Injection properties per Lagrangian time step

            //- Time at start of injection time step [s]
            scalar timeStep0_;

            //- Number of parcels to introduce []
            label nParcels_;

            //- Volume of parcels to introduce [m^3]
            scalar volume_;


    // Protected member functions

        //- Number of parcels to introduce over the time step
        virtual label nParcelsToInject
        (
            const scalar time0,
            const scalar time1
        ) const = 0;

        //- Volume of parcels to introduce over the time step
        virtual scalar volumeToInject
        (
            const scalar time0,
            const scalar time1
        ) const = 0;


public:

    //- Runtime type information
    TypeName("InjectionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        InjectionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct from dictionary
        InjectionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    // Destructor

        virtual ~InjectionModel();


    // Selector

        static autoPtr<InjectionModel<CloudType> > New
        (
            const dictionary& dict,
            CloudType& owner
        );


    // Access

        //- Return const access the owner cloud object
        const CloudType& owner() const;

        //- Return non-const access the owner cloud object for manipulation
        CloudType& owner();

        //- Return the dictionary
        const dictionary& dict() const;

        //- Return the coefficients dictionary
        const dictionary& coeffDict() const;


    // Member Functions

        //- Flag to indicate whether model activates injection model
        virtual bool active() const = 0;


        // Global information

            //- Return the start-of-injection time
            scalar timeStart() const;

            //- Return the total volume to be injected across the event
            scalar volumeTotal() const;

            //- Return the end-of-injection time
            virtual scalar timeEnd() const = 0;


        // Per Lagrangian time step properties

            //- Determine properties for next time step/injection interval
            void prepareForNextTimeStep
            (
                const scalar time0,
                const scalar time1
            );

            //- Return the number of parcels to introduce
            label nParcels() const;

            //- Return the volume of parcels to introduce
            scalar volume() const;

            //- Return the volume fraction to introduce
            scalar volumeFraction() const;


        // Injection geometry

            //- Return the injection position
            virtual vector position
            (
                const label iParcel,
                const scalar time,
                const polyMeshInfo& meshInfo
            ) = 0;

            //- Return the velocity of the parcel to introduce at a time
            virtual vector velocity
            (
                const label iParcel,
                const scalar time,
                const polyMeshInfo& meshInfo
            ) = 0;

            //- Return the diameter of the parcel to introduce at a time
            virtual scalar d0
            (
                const label iParcel,
                const scalar time
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeInjectionModel(CloudType)                                         \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(InjectionModel<CloudType>, 0);        \
                                                                              \
    defineTemplateRunTimeSelectionTable(InjectionModel<CloudType>, dictionary);


#define makeInjectionModelType(SS, CloudType, ParcelType)                     \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType<ParcelType> >, 0);       \
                                                                              \
    InjectionModel<CloudType<ParcelType> >::                                  \
        adddictionaryConstructorToTable<SS<CloudType<ParcelType> > >          \
            add##SS##CloudType##ParcelType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "InjectionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
