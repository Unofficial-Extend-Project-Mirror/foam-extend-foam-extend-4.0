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
    dynamicTopoFvMesh

Description
    An implementation of dynamic changes to mesh-topology

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    dynamicTopoFvMesh.C
    dynamicTopoFvMeshI.H
    dynamicTopoFvMeshCheck.C
    dynamicTopoFvMeshCoupled.C
    dynamicTopoFvMeshReOrder.C
    dynamicTopoFvMeshMapping.C
    edgeBisect.C
    edgeCollapse.C
    edgeSwap.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicTopoFvMesh_H
#define dynamicTopoFvMesh_H

#include "Switch.H"
#include "tetMetric.H"
#include "topoMapper.H"
#include "DynamicField.H"
#include "threadHandler.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class eMesh;
class Stack;
class changeMap;
class objectMap;
class motionSolver;
class convexSetAlgorithm;
class lengthScaleEstimator;
class subMeshLduAddressing;

template <class MeshType>
class coupledInfo;

/*---------------------------------------------------------------------------*\
                      Class dynamicTopoFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicTopoFvMesh
:
    public dynamicFvMesh
{
    // Private typedefs

        typedef coupledInfo<dynamicTopoFvMesh> coupledMesh;
        typedef threadHandler<dynamicTopoFvMesh> meshHandler;

    // Private enumerants

        enum opType
        {
            TOTAL_MODIFICATIONS = 0,
            TOTAL_SWAPS         = 1,
            TOTAL_BISECTIONS    = 2,
            TOTAL_COLLAPSES     = 3,
            SURFACE_SWAPS       = 4,
            SURFACE_BISECTIONS  = 5,
            SURFACE_COLLAPSES   = 6,
            TOTAL_SLIVERS       = 7,
            TOTAL_OP_TYPES
        };

    // Private data

        //- Reference to base mesh
        const dynamicTopoFvMesh& baseMesh_;

        //- Number of boundary patches
        label nPatches_;

        //- Topology change flag
        bool topoChangeFlag_;

        //- Is this a subMesh?
        Switch isSubMesh_;

        //- Dynamic mesh dictionary
        IOdictionary dict_;

        //- Should all options be mandatorily specified?
        //  Handy for first-time use.
        Switch mandatory_;

        //- Mesh characteristics [2D/3D]
        bool twoDMesh_;

        //- Edge refinement switch
        Switch edgeRefinement_;

        //- Load motion solver
        Switch loadMotionSolver_;

        //- Switch for cell-bandwidth reduction
        Switch bandWidthReduction_;

        //- Coupled modification switch
        mutable Switch coupledModification_;

        //- Sub-Mesh lduAddressing
        mutable subMeshLduAddressing* lduPtr_;

        //- Specify the re-meshing interval
        label interval_;

        //- Edge-mesh
        autoPtr<eMesh> eMeshPtr_;

        //- Field mapper
        autoPtr<topoMapper> mapper_;

        //- Mesh motion solver
        autoPtr<motionSolver> motionSolver_;

        //- Length scale estimator
        autoPtr<lengthScaleEstimator> lengthEstimator_;

        //- Lists that dynamically resize during topo-changes
        //   - Since resizes happen infrequently,
        //     scale up by 10% to save memory.
        //   - Current C++ standard doesn't
        //     support template typedefs,
        //     so this is a work-around
        template <class T>
        class resizable
        {
        public:
            typedef DynamicList<T, 0, 11, 10> ListType;
            typedef DynamicField<T, 0, 11, 10> FieldType;
        };

        resizable<point>::FieldType oldPoints_, points_;
        resizable<face>::ListType faces_;
        resizable<label>::ListType owner_, neighbour_;
        resizable<cell>::ListType cells_;
        resizable<edge>::ListType edges_;
        resizable<labelList>::ListType pointEdges_;
        resizable<labelList>::ListType edgeFaces_, faceEdges_;
        resizable<scalar>::ListType lengthScale_;

        //- Size information
        labelList oldPatchSizes_, patchSizes_;
        labelList oldPatchStarts_, patchStarts_;
        labelList oldEdgePatchSizes_, edgePatchSizes_;
        labelList oldEdgePatchStarts_, edgePatchStarts_;
        labelList oldPatchNMeshPoints_, patchNMeshPoints_;
        label nOldPoints_, nPoints_;
        label nOldEdges_, nEdges_;
        label nOldFaces_, nFaces_;
        label nOldCells_, nCells_;
        label nOldInternalFaces_, nInternalFaces_;
        label nOldInternalEdges_, nInternalEdges_;

        //- Maps [Old-to-New]
        labelList reversePointMap_;
        labelList reverseEdgeMap_;
        labelList reverseFaceMap_;
        labelList reverseCellMap_;

        //- Maps [New-to-Old]
        labelList pointMap_;
        labelList edgeMap_;
        labelList faceMap_;
        labelList cellMap_;

        //- Maps for the renumbering of added entities
        Map<label> addedPointRenumbering_;
        Map<label> addedEdgeRenumbering_;
        Map<label> addedFaceRenumbering_;
        Map<label> addedCellRenumbering_;
        Map<label> addedFacePatches_;
        Map<label> addedEdgePatches_;
        Map<label> addedPointZones_;
        Map<label> addedFaceZones_;
        Map<label> addedCellZones_;

        // Information for field-mapping
        Map<labelList> faceParents_;
        List<scalarField> faceWeights_;
        List<vectorField> faceCentres_;

        Map<labelList> cellParents_;
        List<scalarField> cellWeights_;
        List<vectorField> cellCentres_;

        // Information for mapPolyMesh
        List<objectMap> pointsFromPoints_;
        List<objectMap> facesFromPoints_;
        List<objectMap> facesFromEdges_;
        List<objectMap> facesFromFaces_;
        List<objectMap> cellsFromPoints_;
        List<objectMap> cellsFromEdges_;
        List<objectMap> cellsFromFaces_;
        List<objectMap> cellsFromCells_;

        //- Maps to keep track of entities deleted after addition
        labelHashSet deletedPoints_;
        labelHashSet deletedEdges_;
        labelHashSet deletedFaces_;
        labelHashSet deletedCells_;

        //- List of flipped faces
        labelHashSet flipFaces_;

        //- Run-time statistics
        label maxModifications_;
        labelList statistics_;

        //- Sliver exudation
        scalar sliverThreshold_;
        Map<scalar> thresholdSlivers_;

        //- Specific to proximity-based refinement
        List<labelPair> slicePairs_;

        //- 3D Mesh Flipping data
        label maxTetsPerEdge_;
        scalar swapDeviation_;
        Switch allowTableResize_;
        labelList noSwapPatchIDs_;

        //- Stack-list of entities to be checked for topo-changes.
        List<Stack> entityStack_;

        //- Support for multithreading
        autoPtr<IOmultiThreader> threader_;

        // Pointer-list of thread-handlers
        PtrList<meshHandler> handlerPtr_;

        // Entity mutexes used for synchronization
        // in multi-threaded reOrdering
        FixedList<Mutex, 4> entityMutex_;

        // Local coupled patch information
        PtrList<coupledMesh> patchCoupling_;

        // List of processors that this sub-domain talks to
        labelList procIndices_;

        // Processor priority list
        labelList procPriority_;

        // Sub-Mesh pointers
        PtrList<coupledMesh> sendMeshes_;
        PtrList<coupledMesh> recvMeshes_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicTopoFvMesh(const dynamicTopoFvMesh&);

        //- Construct from components
        //- Used for subMeshes only.
        dynamicTopoFvMesh
        (
            const dynamicTopoFvMesh& mesh,
            const IOobject& io,
            const Xfer<pointField>& points,
            const Xfer<pointField>& oldPoints,
            const Xfer<edgeList>& edges,
            const Xfer<faceList>& faces,
            const Xfer<labelListList>& faceEdges,
            const Xfer<labelList>& owner,
            const Xfer<labelList>& neighbour,
            const labelList& faceStarts,
            const labelList& faceSizes,
            const labelList& edgeStarts,
            const labelList& edgeSizes,
            const wordList& patchNames,
            const dictionary& patchDict
        );

        //- Disallow default bitwise assignment
        void operator=(const dynamicTopoFvMesh&);

        // Perform a Delaunay test on an internal face
        bool testDelaunay(const label fIndex) const;

        //- Quality metric for tetrahedra in 3D
        tetMetric::tetMetricReturnType tetMetric_;

        // Compute mapping weights for modified entities
        void computeMapping
        (
            const scalar matchTol,
            const bool skipMapping,
            const bool mappingOutput,
            const label faceStart,
            const label faceSize,
            const label cellStart,
            const label cellSize
        );

        // Static equivalent for multiThreading
        static void computeMappingThread(void *argument);

        // Routine to invoke threaded mapping
        void threadedMapping
        (
            scalar matchTol,
            bool skipMapping,
            bool mappingOutput
        );

        // Initialize mesh edges and related connectivity lists
        void initEdges();

        // Load the mesh-quality metric from the library
        void loadMetric();

        // Load the mesh-motion solver
        void loadMotionSolver();

        // Load the field mapper
        void loadFieldMapper();

        // Load the length scale estimator
        void loadLengthScaleEstimator();

        // Initialize the threading environment
        void initializeThreadingEnvironment(const label specThreads = -1);

        // Return if mesh is 2D
        inline bool is2D() const;

        // Return if mesh is 3D
        inline bool is3D() const;

        // Return a non-const reference to the lengthScaleEstimator
        inline lengthScaleEstimator& lengthEstimator();

        // Return a const reference to the lengthScaleEstimator
        inline const lengthScaleEstimator& lengthEstimator() const;

        // Return a const reference to the multiThreader
        inline const multiThreader& threader() const;

        // Return a reference to the entity mutexes
        inline const Mutex& entityMutex(const label entity) const;

        // Return the edge index for a provided edge
        inline label getEdgeIndex(const edge& edgeToCheck) const;

        // Given a boundary face, pick out a boundary edge that
        // contains a triangular face.
        inline label getTriBoundaryEdge(const label fIndex) const;

        // Given an input quad-face, determine checkEdges from mesh
        static void getCheckEdges
        (
            const label fIndex,
            const dynamicTopoFvMesh& mesh,
            changeMap& map,
            FixedList<edge,4>& checkEdge,
            FixedList<label,4>& checkEdgeIndex
        );

        // Return length-scale at an face-location in the mesh [2D]
        scalar faceLengthScale(const label fIndex) const;

        // Return length-scale at an edge-location in the mesh [3D]
        scalar edgeLengthScale(const label eIndex) const;

        // Check for bisection
        inline bool checkBisection(const label index) const;

        // Check for collapse
        inline bool checkCollapse(const label index) const;

        // MultiThreaded topology modifier
        void threadedTopoModifier();

        // 2D Edge-swapping engine
        static void swap2DEdges(void *argument);

        // 3D Edge-swapping engine
        static void swap3DEdges(void *argument);

        // Edge refinement engine
        static void edgeRefinementEngine(void *argument);

        // Return the entity stack for a particular thread
        inline Stack& stack(const label threadID);

        // Return the integer ID for a given thread
        inline label self() const;

        // Initialize stacks
        inline void initStacks(const labelHashSet& entities);

        // Initialize the coupled stack
        void initCoupledStack
        (
            const labelHashSet& entities,
            bool useEntities
        );

        // Method to determine the boundary patch index for a given face
        inline label whichPatch(const label index) const;

        // Method to determine the boundary patch index for a given edge
        inline label whichEdgePatch(const label index) const;

        // Report or alter topo-modification status
        inline label& status(const label type);

        // Method for the swapping of a quad-face in 2D
        const changeMap
        swapQuadFace
        (
            const label fIndex
        );

        // Method for the bisection of a quad-face in 2D
        const changeMap
        bisectQuadFace
        (
            const label fIndex,
            const changeMap& masterMap,
            bool checkOnly = false,
            bool forceOp = false
        );

        // Method for the collapse of a quad-face in 2D
        const changeMap
        collapseQuadFace
        (
            const label fIndex,
            label overRideCase = -1,
            bool checkOnly = false,
            bool forceOp = false
        );

        // Method for the bisection of an edge in 3D
        const changeMap
        bisectEdge
        (
            const label eIndex,
            bool checkOnly = false,
            bool forceOp = false
        );

        // Method for the collapse of an edge in 3D
        const changeMap
        collapseEdge
        (
            const label eIndex,
            label overRideCase = -1,
            bool checkOnly = false,
            bool forceOp = false
        );

        // Slice the mesh at a particular location
        void sliceMesh(const labelPair& pointPair);

        // Utility method to check for invalid face-bisection.
        bool checkBisection
        (
            const label fIndex,
            const label bFaceIndex,
            bool forceOp = false
        ) const;

        // Utility method to check for invalid face-collapse.
        static bool checkCollapse
        (
            const dynamicTopoFvMesh& mesh,
            const labelList& triFaces,
            const FixedList<label,2>& c0BdyIndex,
            const FixedList<label,2>& c1BdyIndex,
            const FixedList<label,2>& pointIndex,
            const FixedList<point,2>& newPoint,
            const FixedList<point,2>& oldPoint,
            scalar& collapseQuality,
            const bool checkNeighbour,
            bool forceOp = false
        );

        // Utility method to check for invalid edge-collapse.
        bool checkCollapse
        (
            const point& newPoint,
            const point& oldPoint,
            const label pointIndex,
            const label cellIndex,
            dynamicLabelList& cellsChecked,
            scalar& collapseQuality,
            bool forceOp = false
        ) const;

        // Remove 2D sliver cells from the mesh
        void remove2DSlivers();

        // Identify the sliver type in 3D
        const changeMap identifySliverType(const label cIndex) const;

        // Remove sliver cells
        void removeSlivers();

        // Insert the specified cell to the mesh
        label insertCell
        (
            const cell& newCell,
            const scalar lengthScale,
            const label zoneID = -1
        );

        // Set fill-in mapping information for a particular cell
        void setCellMapping
        (
            const label cIndex,
            const labelList& mapCells,
            bool addEntry = true
        );

        // Set fill-in mapping information for a particular face
        void setFaceMapping
        (
            const label fIndex,
            const labelList& mapFaces = labelList(0)
        );

        // Remove the specified cells from the mesh
        const changeMap
        removeCells
        (
            const labelList& cList,
            const label patch,
            const word& rcName,
            bool checkOnly = false
        );

        // Remove the specified cell from the mesh
        void removeCell
        (
            const label cIndex
        );

        // Insert the specified face to the mesh
        label insertFace
        (
            const label patch,
            const face& newFace,
            const label newOwner,
            const label newNeighbour,
            const labelList& newFaceEdges,
            const label zoneID = -1
        );

        // Remove the specified face from the mesh
        void removeFace
        (
            const label fIndex
        );

        // Split a set of internal faces into boundary faces,
        // and add them to the specified patch.
        void splitInternalFaces
        (
            const label patchIndex,
            const labelList& internalFaces,
            const Map<bool>& cellColors
        );

        // Merge a set of boundary faces into internal
        const changeMap
        mergeBoundaryFaces
        (
            const labelList& mergeFaces
        );

        // Insert the specified edge to the mesh
        label insertEdge
        (
            const label patch,
            const edge& newEdge,
            const labelList& edgeFaces
        );

        // Remove the specified edge from the mesh
        void removeEdge
        (
            const label eIndex
        );

        // Insert the specified point to the mesh
        label insertPoint
        (
            const point& newPoint,
            const point& oldPoint,
            const labelList& mapPoints,
            const label zoneID = -1
        );

        // Remove the specified point from the mesh
        void removePoint
        (
            const label pIndex
        );

        // Utility method to build vertexHull for an edge [3D]
        void buildVertexHull
        (
            const label eIndex,
            labelList& vertexHull,
            const label checkIndex = 0
        ) const;

        // Utility to check whether points of an edge lie on a boundary
        const FixedList<bool,2> checkEdgeBoundary(const label eIndex) const;

        // Set a particular face index as flipped.
        inline void setFlip(const label fIndex);

        // Utility method to compute the minimum quality of a vertex hull
        scalar computeMinQuality
        (
            const label eIndex,
            labelList& hullVertices
        ) const;

        // Compute minQuality given addressing
        scalar computeMinQuality
        (
            const edge& edgeToCheck,
            const labelList& hullVertices,
            const UList<point>& points,
            bool closedRing = false
        ) const;

        // Utility method to compute the quality of a vertex hull
        // around an edge after bisection.
        scalar computeBisectionQuality(const label eIndex) const;

        // Check whether the given edge is on a bounding curve
        bool checkBoundingCurve
        (
            const label eIndex,
            const bool overRidePurityCheck = false,
            label* nProcCurves = NULL
        ) const;

        // Allocate dynamic programming tables
        void initTables
        (
            labelList& m,
            PtrList<scalarListList>& Q,
            PtrList<labelListList>& K,
            PtrList<labelListList>& triangulations,
            const label checkIndex = -1
        ) const;

        // Check triangulation quality for an edge index
        bool checkQuality
        (
            const label eIndex,
            const labelList& m,
            const PtrList<scalarListList>& Q,
            const scalar minQuality,
            const label checkIndex = 0
        ) const;

        // Utility method to fill the dynamic programming tables
        bool fillTables
        (
            const label eIndex,
            scalar& minQuality,
            labelList& m,
            labelList& hullVertices,
            PtrList<scalarListList>& Q,
            PtrList<labelListList>& K,
            PtrList<labelListList>& triangulations,
            const label checkIndex = 0
        ) const;

        // Fill tables given addressing
        void fillTables
        (
            const edge& edgeToCheck,
            const scalar minQuality,
            const label m,
            const labelList& hullVertices,
            const UList<point>& points,
            scalarListList& Q,
            labelListList& K,
            labelListList& triangulations
        ) const;

        // Print out tables for debugging
        void printTables
        (
            const labelList& m,
            const PtrList<scalarListList>& Q,
            const PtrList<labelListList>& K,
            const label checkIndex = 0
        ) const;

        // Remove the edge according to the swap sequence
        const changeMap
        removeEdgeFlips
        (
            const label eIndex,
            const scalar minQuality,
            const labelList& hullVertices,
            PtrList<scalarListList>& Q,
            PtrList<labelListList>& K,
            PtrList<labelListList>& triangulations,
            const label checkIndex = 0
        );

        // Extract triangulations from the programming table
        void extractTriangulation
        (
            const label i,
            const label j,
            const labelListList& K,
            label& numTriangulations,
            labelListList& triangulations
        ) const;

        // Identify the 3-2 swap from the triangulation sequence
        label identify32Swap
        (
            const label eIndex,
            const labelList& hullVertices,
            const labelListList& triangulations,
            bool output = false
        ) const;

        // Check old-volumes for the configuration.
        bool checkTriangulationVolumes
        (
            const label eIndex,
            const labelList& hullVertices,
            const labelListList& triangulations
        ) const;

        // Routine to check whether the triangulation at the
        // index lies on the boundary of the vertex ring.
        bool boundaryTriangulation
        (
            const label index,
            label& isolatedVertex,
            labelListList& triangulations
        ) const;

        // Routine to perform 2-3 swaps
        const changeMap
        swap23
        (
            const label isolatedVertex,
            const label eIndex,
            const label triangulationIndex,
            const label numTriangulations,
            const labelListList& triangulations,
            const labelList& hullVertices,
            const labelList& hullFaces,
            const labelList& hullCells
        );

        // Routine to perform 3-2 or 2-2 swaps
        const changeMap
        swap32
        (
            const label eIndex,
            const label triangulationIndex,
            const label numTriangulations,
            const labelListList& triangulations,
            const labelList& hullVertices,
            const labelList& hullFaces,
            const labelList& hullCells
        );

        // Reorder points after a topology change
        void reOrderPoints
        (
            pointField& points,
            pointField& preMotionPoints,
            labelListList& pointZoneMap,
            bool threaded = false
        );

        // Static equivalent for multi-threading
        static void reOrderPointsThread(void *argument);

        // Reorder edges after a topology change
        void reOrderEdges
        (
            edgeList& edges,
            labelListList& edgeFaces,
            labelListList& faceEdges,
            bool threaded = false
        );

        // Static equivalent for multi-threading
        static void reOrderEdgesThread(void *argument);

        // Reorder faces in upper-triangular order after a topology change
        void reOrderFaces
        (
            faceList& faces,
            labelList& owner,
            labelList& neighbour,
            labelListList& faceEdges,
            labelListList& faceZoneFaceMap,
            bool threaded = false
        );

        // Static equivalent for multi-threading
        static void reOrderFacesThread(void *argument);

        // Reorder & renumber cells with bandwidth
        // reduction after a topology change
        void reOrderCells
        (
            labelListList& cellZoneMap,
            bool threaded = false
        );

        // Static equivalent for multi-threading
        static void reOrderCellsThread(void *argument);

        // Reorder the mesh in upper-triangular order,
        // and generate mapping information
        void reOrderMesh
        (
            pointField& points,
            pointField& preMotionPoints,
            edgeList& edges,
            faceList& faces,
            labelList& owner,
            labelList& neighbour,
            labelListList& faceEdges,
            labelListList& edgeFaces,
            labelListList& pointZoneMap,
            labelListList& faceZoneFaceMap,
            labelListList& cellZoneMap
        );

        // Invoke reOrdering with multiple threads
        void threadedMeshReOrdering
        (
            pointField& points,
            pointField& preMotionPoints,
            edgeList& edges,
            faceList& faces,
            labelList& owner,
            labelList& neighbour,
            labelListList& faceEdges,
            labelListList& edgeFaces,
            labelListList& pointZoneMap,
            labelListList& faceZoneFaceMap,
            labelListList& cellZoneMap
        );

        // Reset the mesh and generate mapping information
        bool resetMesh();

        // Write out connectivity for an edge
        void writeEdgeConnectivity(const label eIndex) const;

        // Output an entity as a VTK file
        void writeVTK
        (
            const word& name,
            const label entity,
            const label primitiveType = 3,
            const bool useOldConnectivity = false,
            const bool useOldPoints = false
        ) const;

        // Output a list of entities as a VTK file
        void writeVTK
        (
            const word& name,
            const labelList& cList,
            const label primitiveType = 3,
            const bool useOldConnectivity = false,
            const bool useOldPoints = false,
            const UList<scalar>& scalField = UList<scalar>(),
            const UList<label>& lablField = UList<label>()
        ) const;

        // Return the status report interval
        scalar reportInterval() const;

        // Check the state of local connectivity lists
        void checkConnectivity(const label maxErrors = 0) const;

        // Handle mesh slicing events.
        void handleMeshSlicing();

        // Handle layer addition / removal events
        void handleLayerAdditionRemoval();

        // Add cell layer above specified patch
        const changeMap addCellLayer(const label patchID);

        // Remove cell layer above specified patch
        const changeMap removeCellLayer(const label patchID);

        // Test an edge / face for proximity with other non-neighbouring faces.
        // Return the scalar distance to the nearest face.
        scalar testProximity
        (
            const label index,
            label& proximityFace
        ) const;

        // Calculate the edge length-scale for the mesh
        void calculateLengthScale(bool dump = false);

        // Fill buffers with length-scale info
        // and exchange across processors.
        void exchangeLengthBuffers();

        // Implementing the fillTables operation for coupled edges
        bool coupledFillTables
        (
            const label eIndex,
            scalar& minQuality,
            labelList& m,
            PtrList<scalarListList>& Q,
            PtrList<labelListList>& K,
            PtrList<labelListList>& triangulations
        ) const;

        // Initialize coupled weights calculation
        void initCoupledWeights();

        // Additional mapping contributions for coupled entities
        void computeCoupledWeights
        (
            const label index,
            const label dimension,
            labelList& parents,
            scalarField& weights,
            vectorField& centres,
            bool output = false
        );

        // Fetch length-scale info for processor entities
        scalar processorLengthScale(const label index) const;

        // Method to determine whether the master entity is locally coupled
        bool locallyCoupledEntity
        (
            const label index,
            bool checkSlaves = false,
            bool checkProcs = true,
            bool checkFace = false
        ) const;

        // Method to determine the locally coupled patch index
        label locallyCoupledEdgePatch(const label eIndex) const;

        // Method to determine if the entity is on a processor boundary
        bool processorCoupledEntity
        (
            const label index,
            bool checkFace = false,
            bool checkEdge = false,
            bool checkPure = false,
            FixedList<label, 2>* patchLabels = NULL,
            FixedList<vector, 2>* patchNormals = NULL
        ) const;

        // Read optional dictionary parameters
        void readOptionalParameters(bool reRead = false);

        // Execute load balancing operations on the mesh
        void executeLoadBalancing(dictionary& balanceDict);

        // Identify coupled patches.
        bool identifyCoupledPatches();

        // Read coupled patch information from dictionary.
        void readCoupledPatches();

        // Initialize coupled patch connectivity for topology modifications.
        static void initCoupledConnectivity(void *argument);

        // Move coupled subMesh points
        void moveCoupledSubMeshes();

        // Build a list of entities that need to be avoided
        // by regular topo-changes.
        void buildEntitiesToAvoid
        (
            labelHashSet& entities,
            bool checkSubMesh
        );

        // Check whether the specified edge is a coupled master edge.
        bool isCoupledMaster(const label eIndex) const;

        // Set coupled modification
        void setCoupledModification() const;

        // Unset coupled modification
        void unsetCoupledModification() const;

        // Insert the cells around the coupled master entity to the mesh
        const changeMap insertCells(const label mIndex);

        // Handle topology changes for coupled patches
        void handleCoupledPatches(labelHashSet& entities);

        // Synchronize topology operations across processors
        void syncCoupledPatches(labelHashSet& entities);

        // Check the state of coupled boundaries
        bool checkCoupledBoundaries(bool report = true) const;

        // Build patch sub-meshes for processor patches
        void buildProcessorPatchMeshes();

        // Build patch sub-mesh for a specified processor
        void buildProcessorPatchMesh
        (
            coupledMesh& subMesh,
            Map<labelList>& commonCells
        );

        // Build coupled maps for locally coupled patches
        void buildLocalCoupledMaps();

        // Build coupled maps for coupled processor patches
        void buildProcessorCoupledMaps();

        // Introduce a new processor patch to the mesh
        label createProcessorPatch(const label proc);

        // If a patch is of processor type,
        // get the neighbouring processor ID
        label getNeighbourProcessor(const label patch) const;

        // If the number of patches have changed during run-time,
        // reset boundaries with new processor patches
        void resetBoundaries();

        // Initialize subMesh field transfers for mapping
        void initFieldTransfers
        (
            wordList& fieldTypes,
            List<wordList>& fieldNames,
            List<List<char> >& sendBuffer,
            List<List<char> >& recvBuffer
        );

        // Synchronize field transfers for mapping
        void syncFieldTransfers
        (
            wordList& fieldTypes,
            List<wordList>& fieldNames,
            List<List<char> >& recvBuffer
        );

        // Initialize coupled boundary ordering
        void initCoupledBoundaryOrdering
        (
            List<pointField>& centres,
            List<pointField>& anchors
        ) const;

        // Synchronize coupled boundary ordering
        bool syncCoupledBoundaryOrdering
        (
            List<pointField>& centres,
            List<pointField>& anchors,
            labelListList& faceMaps,
            labelListList& rotations
        ) const;

        // Set processor rank priority
        void initProcessorPriority();

        // Check for processor priority
        template <class BinaryOp>
        inline bool priority
        (
            const label first,
            const BinaryOp& bop,
            const label second
        ) const;

        // Dump cell-quality statistics
        bool meshQuality(bool outputOption);

public:

    //- Runtime type information
    TypeName("dynamicTopoFvMesh");

    // Constructors

        //- Construct from existing IOobject
        explicit dynamicTopoFvMesh(const IOobject& io);

    // Destructor

        virtual ~dynamicTopoFvMesh();

    // Member Functions

        // Return custom lduAddressing
        virtual const lduAddressing& lduAddr() const;

        // Update mesh corresponding to the given map
        virtual void updateMesh(const mapPolyMesh& mpm);

        // Map all fields in time using a customized mapper
        virtual void mapFields(const mapPolyMesh& mpm) const;

        // Update the mesh for motion / topology changes
        //  - Return true if topology changes have occurred
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dynamicTopoFvMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
