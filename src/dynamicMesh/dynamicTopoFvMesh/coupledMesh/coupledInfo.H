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
    coupledInfo

Description
    An interface class that provides patch coupling functionality

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    coupledInfo.C

\*---------------------------------------------------------------------------*/

#ifndef coupledInfo_H
#define coupledInfo_H

#include "autoPtr.H"
#include "coupleMap.H"
#include "volFieldsFwd.H"
#include "fvPatchFieldMapper.H"

namespace Foam
{

// Class forward declarations
class fvBoundaryMesh;

/*---------------------------------------------------------------------------*\
                         Class coupledInfo Declaration
\*---------------------------------------------------------------------------*/

template <class MeshType>
class coupledInfo
{

    // Private data

        // Reference to the parent mesh
        const MeshType& mesh_;

        // Auto pointer to a subMesh
        autoPtr<MeshType> subMesh_;

        // Flag to determine whether maps have been built.
        bool builtMaps_;

        // For locally coupled patches,
        // specify a master / slave index
        coupleMap map_;

        // Zone IDs for patches associated with faceZones
        label masterFaceZone_;
        label slaveFaceZone_;

    // Private member functions

        // Disallow default bitwise assignment
        inline void operator=(const coupledInfo&);

public:

    //- Constructors

        // Construct given mesh, coupleMap and master / slave indices
        inline coupledInfo
        (
            const MeshType& mesh,
            const coupleMap& cMap,
            const label mfzIndex = -1,
            const label sfzIndex = -1
        );

        // Construct from components
        inline coupledInfo
        (
            const MeshType& mesh,
            const bool isTwoDMesh,
            const bool isLocal,
            const bool isSend,
            const label patchIndex,
            const label mPatch,
            const label sPatch,
            const label mfzIndex = -1,
            const label sfzIndex = -1
        );

    //- Access

        // Return a const reference to the parent mesh
        inline const MeshType& baseMesh() const;

        // Set a new subMesh
        inline void setMesh(label index, MeshType* mesh);

        // Return a reference to the subMesh
        inline MeshType& subMesh();

        // Return a const reference to the subMesh
        inline const MeshType& subMesh() const;

        // Return if maps have been built
        inline bool builtMaps() const;

        // Set internal state of maps as built
        inline void setBuiltMaps();

        // Return a reference to the coupleMap
        inline coupleMap& map();

        // Return a const reference to the coupleMap
        inline const coupleMap& map() const;

        // Return the master face zone ID
        inline label masterFaceZone() const;

        // Return the slave face zone ID
        inline label slaveFaceZone() const;

    //- Interpolation

        //- Generic subMesh mapper
        class subMeshMapper
        :
            public fvPatchFieldMapper
        {
            label sizeBeforeMapping_;

            labelField directAddressing_;

        public:

            // Constructors

                //- Construct from components
                inline subMeshMapper
                (
                    const label sbm,
                    const labelList& da
                )
                :
                    sizeBeforeMapping_(sbm),
                    directAddressing_(da)
                {}

                //- Construct given addressing
                inline subMeshMapper
                (
                    const coupledInfo& cInfo,
                    const label patchI
                );

            // Destructor

                virtual ~subMeshMapper()
                {}

            // Member Functions

                label size() const
                {
                    return directAddressing_.size();
                }

                label sizeBeforeMapping() const
                {
                    return sizeBeforeMapping_;
                }

                bool direct() const
                {
                    return true;
                }

                const unallocLabelList& directAddressing() const
                {
                    return directAddressing_;
                }
        };

        // Subset geometric field
        template <class GeomField, class ZeroType>
        tmp<GeomField>
        subSetField
        (
            const GeomField& f,
            const ZeroType& zeroValue,
            const labelList& internalMapper
        ) const;

        // Subset geometric fields from registry to output stream
        template <class GeomField, class ZeroType>
        void send
        (
            const wordList& fieldNames,
            const word& fieldType,
            const ZeroType& zeroValue,
            const labelList& internalMapper,
            OSstream& strStream
        ) const;

        // Set geometric field pointers from input dictionary
        template <class GeomField>
        void setField
        (
            const wordList& fieldNames,
            const dictionary& fieldDicts,
            const label internalSize,
            PtrList<GeomField>& fields
        ) const;

        // Resize map for individual field
        template <class GeomField>
        static void resizeMap
        (
            const label srcIndex,
            const subMeshMapper& internalMapper,
            const List<labelList>& internalReverseMaps,
            const PtrList<subMeshMapper>& boundaryMapper,
            const List<labelListList>& boundaryReverseMaps,
            const List<PtrList<GeomField> >& srcFields,
            GeomField& field
        );

        // Resize map for all fields in registry
        template <class GeomField>
        static void resizeMap
        (
            const wordList& fieldNames,
            const objectRegistry& mesh,
            const subMeshMapper& internalMapper,
            const List<labelList>& internalReverseMaps,
            const PtrList<subMeshMapper>& boundaryMapper,
            const List<labelListList>& boundaryReverseMaps,
            const List<PtrList<GeomField> >& srcFields
        );

        // Resize boundaryFields for all fields in the registry
        template <class GeomField>
        static void resizeBoundaries
        (
            const objectRegistry& mesh,
            const fvBoundaryMesh& boundary
        );
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#    include "coupledInfo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
