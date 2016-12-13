/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::DsmcParcel

Description
    DSMC parcel class

SourceFiles
    DsmcParcelTemplateI.H
    DsmcParcelTemplate.C
    DsmcParcelTemplateIO.C

\*---------------------------------------------------------------------------*/

#ifndef DsmcParcel_H
#define DsmcParcel_H

#include "Particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "contiguous.H"

#include "DsmcCloudTemplate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class DsmcParcel;

// Forward declaration of friend functions

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const DsmcParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                         Class DsmcParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class DsmcParcel
:
    public Particle<ParcelType>
{
public:

    //- Class to hold DSMC particle constant properties
    class constantProperties
    {
        // Private data

            //- Particle mass [kg] (constant)
            scalar mass_;

            //- Particle hard sphere diameter [m] (constant)
            scalar d_;

            //- Internal degrees of freedom
            scalar internalDegreesOfFreedom_;

            //- Viscosity index
            scalar omega_;


    public:

        // Constrcutors

            //- Null constructor, allows List of constantProperties to be
            //  created before the contents is initialised
            inline constantProperties();

            //- Constructor from dictionary
            inline constantProperties(const dictionary& dict);


        // Member functions

            //- Return const access to the particle density
            inline scalar mass() const;

            //- Return const access to the minimum particle mass
            inline scalar d() const;

            //- Return the reference total collision cross section
            inline scalar sigmaT() const;

            //- Return the internalDegreesOfFreedom
            inline scalar internalDegreesOfFreedom() const;

            //- Return the viscosity index
            inline scalar omega() const;

    };


    //- Class used to pass kinematic tracking data to the trackToFace function
    class trackData
    :
        public Particle<ParcelType>::trackData
    {
        // Private data

            //- Reference to the cloud containing this particle
            DsmcCloud<ParcelType>& cloud_;


    public:

        // Constructors

            //- Construct from components
            inline trackData
            (
                DsmcCloud<ParcelType>& cloud
            );


        // Member functions

            //- Return access to the owner cloud
            inline DsmcCloud<ParcelType>& cloud();
    };


protected:

    // Protected member data

        // Parcel properties

            //- Velocity of Parcel [m/s]
            vector U_;

            //- Internal energy of the Parcel, covering all non-translational
            // degrees of freedom [J]
            scalar Ei_;

            //- Parcel type id
            label typeId_;


public:

    //- Runtime type information
    TypeName("DsmcParcel");

    friend class Cloud<ParcelType>;


    // Constructors

        //- Construct from components
        inline DsmcParcel
        (
            DsmcCloud<ParcelType>& owner,
            const vector& position,
            const vector& U,
            const scalar Ei,
            const label celli,
            const label typeId
        );

        //- Construct from Istream
        DsmcParcel
        (
            const Cloud<ParcelType>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<ParcelType> clone() const
        {
            return autoPtr<ParcelType>(new DsmcParcel<ParcelType>(*this));
        }


    // Member Functions

        // Access

            //- Return type id
            inline label typeId() const;

            //- Return const access to velocity
            inline const vector& U() const;

            //- Return const access to internal energy
            inline scalar Ei() const;

        // Edit

            //- Return access to velocity
            inline vector& U();

            //- Return access to internal energy
            inline scalar& Ei();


        // Main calculation loop

            // Tracking

                //- Move the parcel
                template<class TrackData>
                bool move(TrackData& td);


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            template<class TrackData>
            bool hitPatch
            (
                const polyPatch&,
                TrackData& td,
                const label patchI
            );

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            template<class TrackData>
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                TrackData& td
            );

            //- Overridable function to handle the particle hitting a
            //  processorPatch without trackData
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a wallPatch
            template<class TrackData>
            void hitWallPatch
            (
                const wallPolyPatch&,
                TrackData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            //  without trackData
            void hitWallPatch
            (
                const wallPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a polyPatch
            template<class TrackData>
            void hitPatch
            (
                const polyPatch&,
                TrackData& td
            );

            //- Overridable function to handle the particle hitting a polyPatch
            //- without trackData
            void hitPatch
            (
                const polyPatch&,
                int&
            );

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            void transformProperties(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            void transformProperties(const vector& separation);


        // I-O

            static void readFields(Cloud<ParcelType>& c);

            static void writeFields(const Cloud<ParcelType>& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const DsmcParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DsmcParcelTemplateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineParcelTypeNameAndDebug(Type, DebugSwitch)                       \
    template<>                                                                \
    const Foam::word DsmcParcel<Type>::typeName(#Type);                       \
    template<>                                                                \
    Foam::debug::debugSwitch                                                  \
    DsmcParcel<Type>::debug                                                   \
    (                                                                         \
        std::string(#Type), DebugSwitch                                       \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DsmcParcelTemplate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
