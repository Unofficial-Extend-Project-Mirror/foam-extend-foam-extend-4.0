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
    boundaryLayers

Description
    Creates boundary layers

SourceFiles
    boundaryLayers.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryLayers_H
#define boundaryLayers_H

#include "objectRegistry.H"
#include "Time.H"
#include "polyMeshGenModifier.H"
#include "meshSurfaceEngine.H"
#include "meshSurfacePartitioner.H"
#include "DynList.H"
#include "labelLongList.H"
#include "Map.H"
#include "labelPair.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshSurfaceEngine;
class meshSurfacePartitioner;

/*---------------------------------------------------------------------------*\
                           Class boundaryLayers Declaration
\*---------------------------------------------------------------------------*/

class boundaryLayers
{
        //- Reference to the mesh
        polyMeshGen& mesh_;

        //- pointer to mesh surface engine
        mutable meshSurfaceEngine* msePtr_;

        //- poiner to meshSurfacePartitioner
        mutable meshSurfacePartitioner* meshPartitionerPtr_;

        //- shall I create patch-wise layers (true)
        //- O-topology layer (false)
        bool patchWiseLayers_;

        //- shall the layers be terminated at concave edges (true)
        bool terminateLayersAtConcaveEdges_;

        //- is it a 2D mesh
        bool is2DMesh_;

        //- patch names
        wordList patchNames_;

        //- helper which contains information if a boundary layer
        //- has already been extruded for a given patch
        boolList treatedPatch_;

        //- extrude patches with patch
        List<DynList<label> > treatPatchesWithPatch_;

        //- label of a new node (helper)
        labelLongList newLabelForVertex_;

        //- map storing labels of new vertices created at node
        //- and corner vertices
        std::map<label, std::map<std::pair<label, label>, label> > otherVrts_;

        //- a key assigned to each patch. It is needed to search in otherVrts_
        labelList patchKey_;

        //- number of vertices in the mesh
        label nPoints_;

        //- has the geometry been analysed
        bool geometryAnalysed_;

    // Private member functions
        //- Return const reference to meshSurfaceEngine
        const meshSurfaceEngine& surfaceEngine() const;

        //- return const reference to meshSurfacePartitioner
        const meshSurfacePartitioner& surfacePartitioner() const;

        //- find if any other patches need to be treated together
        //- with the given one
        void findPatchesToBeTreatedTogether();

        //- find vertices of the selected patches
        void findPatchVertices
        (
            const boolList& treatPatches,
            List<direction>& patchVertex
        ) const;

        //- check and correct the topology of boundary faces where
        //- the layers terminate
        void checkTopologyOfBoundaryFaces(const labelList& patchLabels);

        //- create new vertex
        point createNewVertex
        (
            const label bpI,
            const boolList& treatPatches,
            const List<direction>& pVertices
        ) const;

        //- create new vertices for the selected patches
        void createNewVertices
        (
            const boolList& treatPatches
        );

        //- create new vertices such that layers for selected patches
        //- are generated in a single run
        void createNewVertices
        (
            const labelList& patchLabels
        );

        //- creates new vertices for vertices at parallel boundaries
        void createNewPartitionVerticesParallel
        (
            const labelLongList& procPoints,
            const List<direction>& pVertices,
            const boolList& treatPatches
        );

        //- creates new vertices from vertices at parallel boundaries
        //- which are also at the border of the treated partitions
        void createNewEdgeVerticesParallel
        (
            const labelLongList& procPoints,
            const List<direction>& pVertices,
            const boolList& treatPatches
        );

        //- create a layer of cells
        void createNewFacesAndCells
        (
            const boolList& treatPatches
        );

        //- create new faces at parallel boundaries
        //- faces are extruded from edges
        void createNewFacesParallel
        (
            const boolList& treatPatches
        );

        //- create new faces at parallel boundaries
        //- faces are created from points at parallel boundaries
        //- the function takes a reference to the faces which are the candidates
        //- to create faces at parallel boundaries
        void createNewFacesFromPointsParallel
        (
            const LongList<DynList<label, 4> >& faceCandidates,
            const LongList<labelPair>& candidatePatches
        );

        //- create layer cells in one go
        //- this is much faster than layer by layer
        void createLayerCells(const labelList& patchLabels);

        //- helper function finding a new face label for multiply extruded nodes
        inline label findNewNodeLabel
        (
            const label pointI,
            const label pKey
        ) const;

        //- creating hex cells near feature edges
        inline void createNewCellFromEdge
        (
            const edge& e,
            const label pKeyI,
            const label pKeyJ,
            FixedList<FixedList<label, 4>, 6>& cellFaces
        ) const;

        //- creating hex cells near corners
        inline void createNewCellFromNode
        (
            const label pointI,
            const DynList<label, 3>& pKeys,
            FixedList<FixedList<label, 4>, 6>& cellFaces
        ) const;

        //- create a bnd layer for a given patch
        void addLayerForPatch(const label patchLabel);

        //- delete meshSurfaceEngine
        inline void clearOut()
        {
            deleteDemandDrivenData(msePtr_);
            deleteDemandDrivenData(meshPartitionerPtr_);
        }

    // Enumerators
        enum vertexTypes
        {
            NONE = 0,
            PATCHNODE = 1,
            EDGENODE = 2,
            CORNERNODE = 4,
            PARALLELBOUNDARY = 8
        };

        //- Disallow bitwise copy construct
        boundaryLayers(const boundaryLayers&);

        //- Disallow bitwise assignment
        void operator=(const boundaryLayers&);

public:

    // Constructors

        //- Construct from mesh reference
        boundaryLayers(polyMeshGen& mesh);

    // Destructor
        ~boundaryLayers();

    // Public member functions
        //- adds layer for a given patch
        void addLayerForPatch(const word& patchName);

        //- create O-topology layers (used as flag)
        void createOTopologyLayers();

        //- terminate boundary layers at concave edges (used as a flag)
        void terminateLayersAtConcaveEdges();

        //- avoid generating layers for empty patches in case of a 2D mesh
        //- used as a flag prior to addLayerForAllPatches
        void activate2DMode();

        //- add layers for all patches
        void addLayerForAllPatches();

        //- add wrapper layer
        //- this function is intended for usage before surface recovery
        void addWrapperLayer();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "boundaryLayersI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
