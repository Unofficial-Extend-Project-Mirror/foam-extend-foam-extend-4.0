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
    coupleMap

Description
    Coupled patch information with registry support.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    coupleMapI.H
    coupleMap.C

\*---------------------------------------------------------------------------*/

#ifndef coupleMap_H
#define coupleMap_H

#include "regIOobject.H"
#include "labelList.H"
#include "pointField.H"
#include "fieldTypes.H"
#include "edgeList.H"
#include "faceList.H"
#include "cellList.H"
#include "labelPair.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class coupleMap Declaration
\*---------------------------------------------------------------------------*/

class coupleMap
:
    public regIOobject
{

public:

    // Public data types

        //- Enumerants for entity-types
        enum entityType
        {
            // Buffers and sizes
            POINT = 0,
            EDGE = 1,
            FACE = 2,
            // Sizes only
            CELL = 3,
            INTERNAL_EDGE = 4,
            INTERNAL_FACE = 5,
            SHARED_POINT = 6,
            GLOBAL_POINT = 7,
            NFE_SIZE = 8,
            NBDY = 9,
            // Buffers only
            OWNER = 3,
            NEIGHBOUR = 4,
            FACE_EDGE = 5,
            NFE_BUFFER = 6,
            FACE_STARTS = 7,
            FACE_SIZES = 8,
            EDGE_STARTS = 9,
            EDGE_SIZES = 10,
            PATCH_ID = 11
        };

        //- Enumerants for operations
        enum opType
        {
            BISECTION = 0,
            COLLAPSE_FIRST = 1,
            COLLAPSE_SECOND = 2,
            COLLAPSE_MIDPOINT = 3,
            REMOVE_CELL = 4,
            MOVE_POINT = 5,
            CONVERT_PATCH = 6,
            INVALID
        };

private:

    // Private data

        mutable bool twoDMesh_;

        // Flags for coupled entities
        mutable bool isLocal_;
        mutable bool isSend_;

        // Patch index for coupled entities
        mutable label patchIndex_;

        // Master / slave indices
        mutable label masterIndex_;
        mutable label slaveIndex_;

        // Point buffers
        mutable pointField pointBuffer_;
        mutable pointField oldPointBuffer_;

        // List of points that are required
        // during patch sub-mesh creation
        mutable labelList subMeshPoints_;
        mutable List<labelPair> globalProcPoints_;

        // Entity sizes (as specified by the entityType enumerant)
        mutable FixedList<label,10> nEntities_;

        // Maps for entities
        mutable FixedList<Map<label>,4> entityMap_;
        mutable FixedList<Map<label>,4> reverseEntityMap_;

        // Entity Buffers (as specified by the entityType enumerant)
        mutable FixedList<labelList,12> entityBuffer_;

        // List of entity indices with topological operations
        mutable labelList entityIndices_;

        // List of operations performed on entities
        mutable List<opType> entityOperations_;

        // List of point-locations to move points to
        mutable pointField moveNewPoints_;
        mutable pointField moveOldPoints_;

        // Addressing for field mapping
        mutable labelList faceMap_;
        mutable labelList cellMap_;

        //- Demand-driven connectivity data.
        mutable edgeList* edgesPtr_;
        mutable faceList* facesPtr_;
        mutable labelListList* faceEdgesPtr_;

        void makeEdges() const;
        void makeFaces() const;

        void makeFaceMap() const;
        void makeCellMap() const;

        void clearAddressing() const;

        //- Disallow default bitwise assignment
        void operator=(const coupleMap&);

public:

    // Static data members

        //- The set of names corresponding to the operationType enumeration
        //  Includes an extra entry for "invalid".
        static const char* names[INVALID + 1];

        //- Runtime type information
        TypeName("coupleMap");

    // Static Member Functions

        //- Return a text representation of an opType
        static const char* asText(const opType);

    // Constructors

        //- Construct from components
        coupleMap
        (
            const IOobject& io,
            const bool twoDMesh,
            const bool isLocal,
            const bool isSend,
            const label patchIndex,
            const label masterIndex,
            const label slaveIndex
        );

        //- Construct as copy
        coupleMap(const coupleMap&);

    // Destructor

        virtual ~coupleMap();

    //- Interpolation functions

        //- Interpolate point field
        template<class Type>
        tmp<Field<Type> > pointInterpolate
        (
            const Map<label>& mPointMap,
            const Map<label>& sPointMap,
            const Field<Type>& pf,
            bool reverse = false
        ) const;

        template<class Type>
        tmp<Field<Type> > pointInterpolate
        (
            const Map<label>& mPointMap,
            const Map<label>& sPointMap,
            const tmp<Field<Type> >& tpf,
            bool reverse = false
        ) const;

        //- Interpolate face field
        template<class Type>
        tmp<Field<Type> > faceInterpolate
        (
            const label mStart,
            const label sStart,
            const Field<Type>& pf,
            bool reverse = false
        ) const;

        template<class Type>
        tmp<Field<Type> > faceInterpolate
        (
            const label mStart,
            const label sStart,
            const tmp<Field<Type> >& tpf,
            bool reverse = false
        ) const;

    //- Access

        inline label& patchIndex() const;

        inline label masterIndex() const;
        inline label slaveIndex() const;

        inline bool isTwoDMesh() const;
        inline bool isLocal() const;
        inline bool isProcessor() const;
        inline bool isSend() const;
        inline bool isRecv() const;

        pointField& pointBuffer() const;
        pointField& oldPointBuffer() const;

        labelList& subMeshPoints() const;
        List<labelPair>& globalProcPoints() const;

        void allocateBuffers() const;

        label findSlave
        (
            const label eType,
            const label Index
        ) const;

        label findMaster
        (
            const label eType,
            const label Index
        ) const;

        void removeSlave
        (
            const label eType,
            const label Index
        ) const;

        void removeMaster
        (
            const label eType,
            const label Index
        ) const;

        void mapSlave
        (
            const label eType,
            const label master,
            const label slave
        ) const;

        void mapMaster
        (
            const label eType,
            const label slave,
            const label master
        ) const;

        inline FixedList<label,10>& nEntities() const;
        inline label& nEntities(const label eType) const;

        inline Map<label>& entityMap(const label eType) const;
        inline Map<label>& reverseEntityMap(const label eType) const;

        inline FixedList<labelList,12>& entityBuffer() const;
        inline labelList& entityBuffer(const label eType) const;

        inline labelList& entityIndices() const;
        inline List<opType>& entityOperations() const;

        inline pointField& moveNewPoints() const;
        inline pointField& moveOldPoints() const;

        void pushOperation
        (
            const label index,
            const opType oType,
            const point& newPoint = vector::zero,
            const point& oldPoint = vector::zero
        ) const;

        //- Demand-driven data
        label nInternalFaces() const;
        const labelList& owner() const;
        const labelList& neighbour() const;
        const edgeList& edges() const;
        const faceList& faces() const;
        const labelListList& faceEdges() const;

        const labelList& faceMap() const;
        const labelList& cellMap() const;

        void transferMaps
        (
            const label eType,
            Map<label>& newEntityMap,
            Map<label>& newReverseEntityMap
        ) const;

        void clearMaps() const;

        void clearBuffers() const;

        bool readData(Istream&);

        bool writeData(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "coupleMapI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
