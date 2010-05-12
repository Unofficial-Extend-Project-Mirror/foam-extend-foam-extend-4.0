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
    Foam::Cloud

Description

SourceFiles
    Cloud.C
    CloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef Cloud_H
#define Cloud_H

#include "cloud.H"
#include "pointMesh.H"
#include "IDLList.H"
#include "IOField.H"

#include "polyMeshInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParticleType>
class Cloud;

template<class ParticleType>
Ostream& operator<<
(
    Ostream&,
    const Cloud<ParticleType>&
);


/*---------------------------------------------------------------------------*\
                           Class Cloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class Cloud
:
    public cloud,
    public IDLList<ParticleType>
{

    // Private data

        const polyMesh& polyMesh_;
        const faceList& allFaces_;
        const vectorField& points_;
        const cellList& cellFaces_;
        const vectorField& allFaceCentres_;
        const unallocLabelList& owner_;
        const unallocLabelList& neighbour_;

        //- Mesh information object
        const polyMeshInfo meshInfo_;


    // Private member functions

        //- Initialise cloud on IO constructor
        void initCloud(const bool checkClass);


public:

    template<class ParticleT>
    friend class Particle;

    typedef ParticleType particleType;

    typedef typename IDLList<ParticleType>::iterator iterator;
    typedef typename IDLList<ParticleType>::const_iterator const_iterator;

    //-Runtime type information
    TypeName("Cloud");


    // Constructors

        //- Construct from mesh and a list of particles
        Cloud
        (
            const polyMesh& mesh,
            const IDLList<ParticleType>& particles
        );

        //- Construct from mesh, cloud name, and a list of particles
        Cloud
        (
            const polyMesh& mesh,
            const word& cloudName,
            const IDLList<ParticleType>& particles
        );

        //- Construct from mesh by reading from file
        //  Optionally disable checking of class name for post-processing
        Cloud
        (
            const polyMesh& mesh,
            const bool checkClass = true
        );


        //- Construct from mesh by reading from file with given cloud instance
        //  Optionally disable checking of class name for post-processing
        Cloud
        (
            const polyMesh& pMesh,
            const word& cloudName,
            const bool checkClass = true
        );


    // Member Functions

        // Access

            //- Return the polyMesh reference
            const polyMesh& pMesh() const
            {
                return polyMesh_;
            }

            //- Is this global face an internal face?
            bool internalFace(const label facei) const
            {
                return polyMesh_.isInternalFace(facei);
            }

            //- Is this global face a boundary face?
            bool boundaryFace(const label facei) const
            {
                return !internalFace(facei);
            }

            //- Which patch is this global face on
            label facePatch(const label facei) const
            {
                return polyMesh_.boundaryMesh().whichPatch(facei);
            }

            //- Which face of this patch is this global face
            label patchFace(const label patchi, const label facei) const
            {
                return polyMesh_.boundaryMesh()[patchi].whichFace(facei);
            }

            //- Return reference to the mesh information object
            const polyMeshInfo& meshInfo() const
            {
                return meshInfo_;
            }

            label size() const
            {
                return IDLList<ParticleType>::size();
            };

            const const_iterator begin() const
            {
                return IDLList<ParticleType>::begin();
            };

            const const_iterator end() const
            {
                return IDLList<ParticleType>::end();
            };

            iterator begin()
            {
                return IDLList<ParticleType>::begin();
            };

            iterator end()
            {
                return IDLList<ParticleType>::end();
            };

            void clear()
            {
                return IDLList<ParticleType>::clear();
            };


        // Edit

            //- Transfer particle to cloud
            void addParticle(ParticleType* pPtr);

            //- Remove particle from cloud and delete
            void deleteParticle(ParticleType&);

            //- Move the particles
            //  passing the TrackingData to the track function
            template<class TrackingData>
            void move(TrackingData& td);

            //- Remap the cells of particles corresponding to the
            //  mesh topology change
            virtual void autoMap(const mapPolyMesh&);


        // Read

            //- Read and return a lagrangian data field
            //  Checks that data is valid
            template<class DataType>
            void checkFieldIOobject
            (
                const Cloud<ParticleType>& c,
                const IOField<DataType>& data
            ) const;

            //- Read and return a lagrangian data field
            IOobject fieldIOobject(const word& fieldName) const;

            //- Read and return a lagrangian data field
            template<class Type>
            tmp<IOField<Type> > readField(const word& fieldName) const;

            //- Read the field data for the cloud of particles
            void readFields();


        // Write

            //- Write the field data for the cloud of particles
            virtual void writeFields() const;

            //- Write using given format, version and compression.
            //  Only writes the cloud file if the Cloud isn't empty
            virtual bool writeObject
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp
            ) const;


    // Ostream Operator

        friend Ostream& operator<< <ParticleType>
        (
            Ostream&,
            const Cloud<ParticleType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Cloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
