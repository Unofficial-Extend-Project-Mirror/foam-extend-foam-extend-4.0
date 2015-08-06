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
    partTriMesh

Description
    Triangulation of mesh surface needed for surface smoothing

SourceFiles
    partTriMesh.C

\*---------------------------------------------------------------------------*/

#ifndef partTriMesh_H
#define partTriMesh_H

#include "boolList.H"
#include "labelLongList.H"
#include "VRWGraph.H"
#include "DynList.H"
#include "HashSet.H"
#include "labelledPoint.H"
#include "triSurf.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshSurfacePartitioner;
class VRWGraph;

/*---------------------------------------------------------------------------*\
                        Class partTriMesh Declaration
\*---------------------------------------------------------------------------*/

class partTriMesh
{
    // Private data
        //- const reference to the meshSurfacePartitioner
        const meshSurfacePartitioner& mPart_;

        //- surface triangulation created from
        triSurf surf_;

        //- label of point in the mesh surface
        labelLongList pointLabelInMeshSurface_;

        //- label of mesh surface point in the partTriMesh
        labelList meshSurfacePointLabelInTriMesh_;

        //- shall a node be used for smoothing or not
        LongList<direction> pointType_;

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

        //- create surface triangulation
        void createPointsAndTrias(const List<direction>&);

        //- create parallel addressing
        void createParallelAddressing
        (
            const labelList& nodeLabelForPoint,
            const labelList& nodeLabelForFace
        );

        //- create buffer layers
        void createBufferLayers();

        //- update buffer layer points
        void updateBufferLayers();

        //- disallow bitwise assignment
        void operator=(const partTriMesh&);

        //- disallow bitwise copy construct
        partTriMesh(const partTriMesh&);

public:

    // Constructors
        //- construct from meshSurfacePartitioner
        partTriMesh(const meshSurfacePartitioner& mPart);

        //- construct from meshSurfacePartitioner, inverted points
        //- and the number of additional layers
        partTriMesh
        (
            const meshSurfacePartitioner& mPart,
            const labelHashSet& invertedPoints,
            const label additionalLayers = 0
        );

    // Enumerators

        enum vertexTypes
        {
            NONE = 0,
            SMOOTH = 1,
            FACECENTRE = 2,
            PARALLELBOUNDARY = 8,
            BOUNDARY = 16,
            FEATUREEDGE = 32,
            CORNER = 64
        };

    // Destructor
        ~partTriMesh();

    // Member functions
        //- access to points, tets and other data
        inline const pointField& points() const
        {
            return surf_.points();
        }

        inline const LongList<labelledTri>& triangles() const
        {
            return surf_.facets();
        }

        inline const VRWGraph& pointTriangles() const
        {
            return surf_.pointFacets();
        }

        inline const LongList<direction>& pointType() const
        {
            return pointType_;
        }

        inline const labelLongList& pointLabelInMeshSurface() const
        {
            return pointLabelInMeshSurface_;
        }

        //- return indices of mesh sutrface points in the surface triangulation
        //- additional points which do not exist in mesh surface are labelled -1
        inline const labelList& meshSurfacePointLabelInTriMesh() const
        {
            return meshSurfacePointLabelInTriMesh_;
        }

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

        //- update coordinates of points in partTriMesh to match the coordinates
        //- in the mesh surface
        void updateVertices();

        //- update coordinates of points in partTriMesh to match the coordinates
        //- of the specified points in the mesh surface
        void updateVertices(const labelLongList&);

        //- return triSurf from this partTriMesh
        const triSurf& getTriSurf() const
        {
            return surf_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
