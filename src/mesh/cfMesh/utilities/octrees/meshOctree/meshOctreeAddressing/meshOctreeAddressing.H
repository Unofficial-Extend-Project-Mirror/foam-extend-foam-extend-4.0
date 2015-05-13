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
    meshOctreeAddressing

Description
    Automatic octree for mesh generation. This class checks if the current
    octree resolution is enough to resolve all surface features

SourceFiles
    meshOctreeAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef meshOctreeAddressing_H
#define meshOctreeAddressing_H

#include "FRWGraph.H"
#include "VRWGraph.H"
#include "meshOctreeCubeCoordinates.H"
#include "face.H"
#include "edge.H"
#include "DynList.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class IOdictionary;

/*---------------------------------------------------------------------------*\
                    Class meshOctreeAddressing Declaration
\*---------------------------------------------------------------------------*/

class meshOctreeAddressing
{
    // Private data
        //- reference to the octree
        const meshOctree& octree_;
    
        //- reference to the dictionary
        const dictionary& meshDict_;

        //- use DATA boxes
        bool useDATABoxes_;
    
        //- number of created octree nodes
        mutable label nNodes_;
    
        //- coordinates of octree nodes
        mutable pointField* octreePointsPtr_;
    
        //- node labels
        mutable VRWGraph* nodeLabelsPtr_;
    
        //- node leaves
        mutable FRWGraph<label, 8>* nodeLeavesPtr_;
    
        //- identify which boxes should be used as mesh cells
        mutable List<direction>* boxTypePtr_;
        
        //- identify created nodes as OUTERNODE or INNERNODE
        mutable List<direction>* nodeTypePtr_;
        
    // Additional addressing useful for some algorithms
        //- faces of the octree
        mutable VRWGraph* octreeFacesPtr_;
        mutable labelLongList* octreeFacesOwnersPtr_;
        mutable labelLongList* octreeFacesNeighboursPtr_;
        
        //- octree box-faces addressing
        mutable VRWGraph* leafFacesPtr_;
        
        //- node-faces addressing
        mutable VRWGraph* nodeFacesPtr_;
        
        //- leaf-leaves addressing
        mutable VRWGraph* leafLeavesPtr_;
        
        //- edges of the octree
        mutable LongList<edge>* octreeEdgesPtr_;
        
        //- edges-leaves addressing
        mutable VRWGraph* edgeLeavesPtr_;
        
        //- leaf-edges addressing
        mutable VRWGraph* leafEdgesPtr_;
        
        //- node-edges addressing
        mutable VRWGraph* nodeEdgesPtr_;
        
        //- face-edges addressing
        mutable VRWGraph* faceEdgesPtr_;
        
        //- edge-faces addressing
        mutable VRWGraph* edgeFacesPtr_;
        
        //- global octree point label
        mutable labelLongList* globalPointLabelPtr_;
        
        //- global point to local label addressing
        mutable Map<label>* globalPointToLocalPtr_;
        
        //- point-processors addressing
        mutable VRWGraph* pointProcsPtr_;
        
        //- global octree face label
        mutable labelLongList* globalFaceLabelPtr_;
        
        //- global face label to local label addressing
        mutable Map<label>* globalFaceToLocalPtr_;
        
        //- face-processors addressing
        mutable VRWGraph* faceProcsPtr_;
        
        //- global leaf label
        mutable labelLongList* globalLeafLabelPtr_;
        
        //- global leaf label to local label addressing for octree leaves
        mutable Map<label>* globalLeafToLocalPtr_;
        
        //- leaf at procs
        mutable VRWGraph* leafAtProcsPtr_;
    
    // Private member functions
        //- calculate octreePointsPtr_
        void createOctreePoints() const;
        
        //- calculate nodeLabelsPtr_
        void createNodeLabels() const;
    
        //- calculate nodeLeavesPtr_
        void createNodeLeaves() const;
        
        //- assemble boxTypePtr_ list
        void findUsedBoxes() const;
        
        //- calculate nodeTypePtr_
        void calculateNodeType() const;
        
        //- calculate faces
        void createOctreeFaces() const;
        
        //- calculate leaf-faces addressing
        void calculateLeafFaces() const;
        
        //- calculate node-faces addressing
        void calculateNodeFaces() const;
        
        //- calculate leaf-leaves addressing
        void calculateLeafLeaves() const;
        
        //- calculate edges
        void createOctreeEdges() const;
        
        //- calculate leaf-edges
        void calculateLeafEdges() const;
        
        //- calculate edge-leaves addressing
        void calculateEdgeLeaves() const;
        
        //- calculate edge-faces
        void calculateEdgeFaces() const;
    
        //- Clear allocated data
        void clearOut();
        void clearBoxTypes();
        void clearNodeAddressing();
        void clearOctreeFaces();
        void clearAddressing();
        void clearParallelAddressing();
        
        //- check if distinct parts are glued together
        void checkGluedRegions();
    
    // Private member functions for parallel runs
    
        void calcGlobalPointLabels() const;
        void calcGlobalFaceLabels() const;
        void calcGlobalLeafLabels() const;

    // Private copy constructor
        //- Disallow default bitwise copy construct
        meshOctreeAddressing(const meshOctreeAddressing&);

        //- Disallow default bitwise assignment
        void operator=(const meshOctreeAddressing&);

public:

    // Constructors

        //- Construct from surface and dictionary
        meshOctreeAddressing
        (
            const meshOctree& mo,
            const dictionary& dict,
            bool useDATABoxes = false
        );

    // Destructor

        ~meshOctreeAddressing();

    // Enumerators

        enum octreeCubeTypes
        {
            NONE = 0,
            MESHCELL = 1,
            BOUNDARY = 2,
            SPLITHEX = 4,
            INNERNODE = 8,
            OUTERNODE = 16,
            BOUNDARYNODE = 32
        };
        
    // Member Functions
        
        //- check if there exist any non-signly connected edges and vertices
        //- in the octree mesh. This is primarily needed
        //- for the cartesian template
        void checkAndFixIrregularConnections();

        //- return number of octree nodes
        inline label numberOfNodes() const;
        
        //- return coordinates of octree vertices
        inline const pointField& octreePoints() const;

        //- return nodeLabels
        inline const VRWGraph& nodeLabels() const;

        //- return nodeLeaves
        inline const FRWGraph<label, 8>& nodeLeaves() const;
        
        //- return which octree boxes are used for mesh creation
        inline const List<direction>& boxType() const;
        
        //- return type of node (INNERNODE,or OUTERNODE)
        inline const List<direction>& nodeType() const;
        
        //- set box type
        inline void setBoxType(const label boxI, const direction type);
        
        //- return octree faces, created for MESHCELL boxes
        inline const VRWGraph& octreeFaces() const;
        
        //- return owners of octree faces
        inline const labelLongList& octreeFaceOwner() const;
        
        //- return neighbours of octree faces
        inline const labelLongList& octreeFaceNeighbour() const;
        
        //- return octree box-faces addressing
        inline const VRWGraph& leafFaces() const;
        
        //- checks if the face is intersected by the surface
        bool isIntersectedFace(const label fI) const;
        
        //- return leaf-leaves addressing
        inline const VRWGraph& leafLeaves() const;
        
        //- return node-faces addressing
        inline const VRWGraph& nodeFaces() const;
        
        //- return octree edges, created for MESHCELL boxes
        inline const LongList<edge>& octreeEdges() const;
        
        //- return edge-leaves addressing
        inline const VRWGraph& edgeLeaves() const;
        
        //- return leaf-edges addressing
        inline const VRWGraph& leafEdges() const;
        
        //- checks if the edge is intersected by the surface
        bool isIntersectedEdge(const label eI) const;
        void edgeIntersections
        (
            const label eI,
            DynList<point>& intersections
        ) const;
        
        //- return node-edges addressing
        inline const VRWGraph& nodeEdges() const;
        
        //- return face-edges addressing
        inline const VRWGraph& faceEdges() const;
        
        //- return edge-faces addressing
        inline const VRWGraph& edgeFaces() const;
        
        //- return const reference to meshOctree
        inline const meshOctree& octree() const;
        
        //- find cubes around an edge (cubes must be at the same level)
        void cubesAroundEdge
        (
            const label leafI,
            const direction eI,
            FixedList<label, 4>& edgeCubes
        ) const;
        
        //- find edge centre if it exists
        label findEdgeCentre(const label leafI, const direction eI) const;
    
    // Access to data needed for parallel execution
    
        //- return global point labels
        inline const labelLongList& globalPointLabel() const;
        
        //- global point label to local label. Only for processors points
        inline const Map<label>& globalToLocalPointAddressing() const;
        
        //- processors which contain an octree point
        inline const VRWGraph& pointAtProcs() const;
        
        //- return global labels of octree faces
        inline const labelLongList& globalFaceLabel() const;
        
        //- return global face label to face label. Only for processor faces
        inline const Map<label>& globalToLocalFaceAddressing() const;
        
        //- return processors which contain each octree face
        inline const VRWGraph& faceAtProcs() const;
        
        //- return global labels of octree leaves
        inline const labelLongList& globalLeafLabel() const;
        
        //- return processors which contain each octree leaf
        inline const VRWGraph& leafAtProcs() const;
        
        //- return global leaf label to local label addressing
        inline const Map<label>& globalToLocalLeafAddressing() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshOctreeAddressingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
