/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Class
    partTetMesh

Description
    Mesh smoothing without any topological changes

SourceFiles
    partTetMesh.C

\*---------------------------------------------------------------------------*/

#ifndef partTetMesh_H
#define partTetMesh_H

#include "boolList.H"
#include "labelLongList.H"
#include "VRWGraph.H"
#include "DynList.H"
#include "partTet.H"
#include "HashSet.H"
#include "labelledPoint.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyMeshGen;
class VRWGraph;

/*---------------------------------------------------------------------------*\
                        Class partTetMesh Declaration
\*---------------------------------------------------------------------------*/

class partTetMesh
{
    // Private data
        //- reference to the original mesh
        polyMeshGen& origMesh_;

        //- points in the tet mesh
        LongList<point> points_;

        //- tetrahedra making the mesh
        LongList<partTet> tets_;

        //- label of node in the polyMeshGen
        labelLongList nodeLabelInOrigMesh_;

        //- shall a node be used for smoothing or not
        LongList<direction> smoothVertex_;

        //- addressing data
        VRWGraph pointTets_;

        //- internal point ordering for parallel runs
        mutable VRWGraph* internalPointsOrderPtr_;

        //- boundary point ordering for parallel runs
        mutable VRWGraph* boundaryPointsOrderPtr_;

    // Private data for parallel runs
        //- global point label
        mutable labelLongList* globalPointLabelPtr_;

        //- processor for containing points
        mutable VRWGraph* pAtProcsPtr_;

        //- mapping between global and local point labels
        mutable Map<label>* globalToLocalPointAddressingPtr_;

        //- processors which should communicate with the current one
        mutable DynList<label>* neiProcsPtr_;

        //- labels of points at parallel boundaries
        mutable labelLongList* pAtParallelBoundariesPtr_;

        //- labels of points serving as buffer layers on other processors
        mutable labelLongList* pAtBufferLayersPtr_;

    // Private member functions

        //- create points and tets
        void createPointsAndTets
        (
            const List<direction>& usedCells,
            const boolList& lockedPoints
        );

        //- create parallel addressing
        void createParallelAddressing
        (
            const labelLongList& nodeLabelForPoint,
            const labelLongList& nodeLabelForFace,
            const labelLongList& nodeLabelForCell
        );

        //- create buffer layers
        void createBufferLayers();

        //- create ordering of internal points for parallel execution
        void createSMOOTHPointsOrdering() const;

        //- create order of boundary points for parallel execution
        void createBOUNDARYPointsOrdering() const;

public:

    // Constructors
        //- construct from polyMeshGen and locked points
        partTetMesh(polyMeshGen& mesh, const labelLongList& lockedPoints);

        //- construct from polyMeshGen, locked points and the number of layers
        //- from the boundary
        partTetMesh
        (
            polyMeshGen& mesh,
            const labelLongList& lockedPoints,
            const direction nLayers = 2
        );

        //- construct from polyMeshGen, bad faces and the number
        //- of additional layers
        partTetMesh
        (
            polyMeshGen& mesh,
            const labelLongList& lockedPoints,
            labelHashSet& badFaces,
            const direction additionalLayers = 0
        );

    // Enumerators

        enum vertexTypes
        {
            NONE = 0,
            SMOOTH = 1,
            FACECENTRE = 2,
            CELLCENTRE = 4,
            PARALLELBOUNDARY = 8,
            BOUNDARY = 16,
            LOCKED = 32
        };

    // Destructor
        ~partTetMesh();

    // Member functions
        //- access to points, tets and other data
        inline const LongList<point>& points() const
        {
            return points_;
        }

        inline const LongList<partTet>& tets() const
        {
            return tets_;
        }

        inline const VRWGraph& pointTets() const
        {
            return pointTets_;
        }

        inline const LongList<direction>& smoothVertex() const
        {
            return smoothVertex_;
        }

        inline const labelLongList& nodeLabelInOrigMesh() const
        {
            return nodeLabelInOrigMesh_;
        }

    // Access to point ordering for Gauss-Seidel type of iteration
        //- return vertex ordering which can be executed in parallel
        //- the points in each row can be treated in parallel
        //- make sure that points added to different rows of this graph are
        //- not treated concurrently
        const VRWGraph& internalPointOrdering() const;

        //- return vertex ordering for boundary points
        const VRWGraph& boundaryPointOrdering() const;

    // Access to parallel data
        inline const labelLongList& globalPointLabel() const
        {
            if( !Pstream::parRun() )
                FatalError << "This is a serial run" << abort(FatalError);

            return *globalPointLabelPtr_;
        }

        inline const VRWGraph& pointAtProcs() const
        {
            if( !Pstream::parRun() )
                FatalError << "This is a serial run" << abort(FatalError);

            return *pAtProcsPtr_;
        }

        inline const Map<label>& globalToLocalPointAddressing() const
        {
            if( !Pstream::parRun() )
                FatalError << "This is a serial run" << abort(FatalError);

            return *globalToLocalPointAddressingPtr_;
        }

        inline const DynList<label>& neiProcs() const
        {
            if( !Pstream::parRun() )
                FatalError << "This is a serial run" << abort(FatalError);

            return *neiProcsPtr_;
        }

        inline const labelLongList& pointsAtProcessorBoundaries() const
        {
            if( !Pstream::parRun() )
                FatalError << "This is a serial run" << abort(FatalError);

            return *pAtParallelBoundariesPtr_;
        }

        inline const labelLongList& bufferLayerPoints() const
        {
            if( !Pstream::parRun() )
                FatalError << "This is a serial run" << abort(FatalError);

            return *pAtBufferLayersPtr_;
        }

    // Modifiers

        //- move the vertex to a new position
        void updateVertex(const label pointI, const point& newP);

        //- move vertices to their new positions
        //- intended for SMP parallelisation
        void updateVerticesSMP(const List<LongList<labelledPoint> >&);

        //- updates the vertices of the original polyMeshGen
        void updateOrigMesh(boolList* changedFacePtr = NULL);

        //- creates polyMeshGen from this partTetMesh
        void createPolyMesh(polyMeshGen& pmg) const;

        //- add the LOCKED flag to mesh points which shall not move
        template<class labelListType>
        void lockPoints(const labelListType& l)
        {
            forAll(l, pI)
                smoothVertex_[l[pI]] |= LOCKED;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
