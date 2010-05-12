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
    Foam::ConeInjection

Description
    Cone injection

    - User specifies
      - time of start of injection
      - injector position
      - direction (along injection axis)
      - parcel flow rate
      - parcel velocity
      - inner and outer cone angles
    - Parcel diameters obtained by PDF model

SourceFiles
    ConeInjection.C

\*---------------------------------------------------------------------------*/

#ifndef ConeInjection_H
#define ConeInjection_H

#include "InjectionModel.H"
#include "pdf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class DataEntry;

/*---------------------------------------------------------------------------*\
                       Class ConeInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ConeInjection
:
    public InjectionModel<CloudType>
{

    // Private data

        //- Coefficients dictionary
        dictionary coeffDict_;

        //- Injection duration [s]
        const scalar duration_;

        //- Injector position [m]
        const vector position_;

        //- Injector direction []
        vector direction_;

        //- Number of parcels to introduce per second []
        const label parcelsPerSecond_;

        //- Volume flow rate of parcels to introduce relative to SOI [m^3]
        const autoPtr<DataEntry<scalar> > volumeFlowRate_;

        //- Parcel velocity magnitude relative to SOI [m/s]
        const autoPtr<DataEntry<scalar> > Umag_;

        //- Inner cone angle relative to SOI [deg]
        const autoPtr<DataEntry<scalar> > thetaInner_;

        //- Outer cone angle relative to SOI [deg]
        const autoPtr<DataEntry<scalar> > thetaOuter_;

        //- Parcel size PDF model
        const autoPtr<pdf> parcelPDF_;


        // Tangential vectors to the direction vector

            //- First tangential vector
            vector tanVec1_;

            //- Second tangential vector
            vector tanVec2_;


protected:

    // Protected member functions

        //- Number of parcels to introduce over the time step
        label nParcelsToInject
        (
            const scalar time0,
            const scalar time1
        ) const;

        //- Number of parcels to introduce over the time step
        scalar volumeToInject
        (
            const scalar time0,
            const scalar time1
        ) const;


public:

    //- Runtime type information
    TypeName("ConeInjection");


    // Constructors

        //- Construct from dictionary
        ConeInjection
        (
            const dictionary& dict,
            CloudType& owner
        );


    // Destructor

        ~ConeInjection();


    // Member Functions

        //- Flag to indicate whether model activates injection model
        bool active() const;

        //- Return the end-of-injection time
        scalar timeEnd() const;


        // Injection geometry

            //- Return the injection position
            vector position
            (
                const label iParcel,
                const scalar time,
                const polyMeshInfo& meshInfo
            );

            //- Return the velocity of the parcel to introduce at a time
            vector velocity
            (
                const label,
                const scalar time,
                const polyMeshInfo& meshInfo
            );

            //- Return the diameter of the parcel to introduce at a time
            scalar d0
            (
                const label,
                const scalar
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ConeInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
