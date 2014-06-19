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

Namespace
    edgeExtractor

Description
    A class with a functions used to detect feature edges in the surface
    of the volume mesh and to detect the patches to which the boundary
    faces belong to

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef edgeExtractor_H
#define edgeExtractor_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "objectRegistry.H"
#include "Time.H"
#include "labelList.H"
#include "labelLongList.H"
#include "VRWGraph.H"
#include "boolList.H"
#include "labelledPoint.H"
#include "DynList.H"
#include "Map.H"

namespace Foam
{

// Forward declarations
class polyMeshGen;
class meshSurfaceEngine;
class meshSurfacePartitioner;
class meshOctree;
class triSurf;
class triSurfaceClassifyEdges;
class triSurfacePartitioner;

/*---------------------------------------------------------------------------*\
                    Class edgeExtractor Declaration
\*---------------------------------------------------------------------------*/

class edgeExtractor
{
    // Private data
        //- reference to the mesh
        polyMeshGen& mesh_;

        //- surface engine
        mutable meshSurfaceEngine* surfaceEnginePtr_;

        //- reference to the octree
        const meshOctree& meshOctree_;

        //- surface mesh partitioner
        mutable const triSurfacePartitioner* surfPartitionerPtr_;

        //- classification of edges in the surface mesh
        mutable const triSurfaceClassifyEdges* surfEdgeClassificationPtr_;

        //- valence of surface points
        labelLongList pointValence_;

        //- patch to which a boundary point is mapped to
        labelLongList pointPatch_;

        //- boundary face patch
        labelList facePatch_;

        //- number of cells attached to a boundary edge
        labelLongList nCellsAtEdge_;

        //- edge classification
        LongList<direction> edgeType_;

        //- patches in the vicinity of a face on the surface of the volume mesh
        VRWGraph patchesNearFace_;

        //- feature edges in the surface mesh
        //- which are candindates for an edge at the surface of the volume mesh
        VRWGraph featureEdgesNearEdge_;

    // Private member functions
        //- calculate point valence
        void calculateValence();

        //- calculate the number of boundary faces for each cell
        void calculateSingleCellEdge();

        //- search for surface patches in the vicinity of a face
        //- at the surface of the volume mesh
        void findPatchesNearSurfaceFace();

        //- search for feature edges in the surface mesh which are
        //- in the vicinity of an edge at the surface of the volume mesh
        void findFeatureEdgesNearEdge();

        //- calculate a list of booleans with the following properties:
        //- true if all faces at a surface vertex are in the same patch
        //- false otherwise
        void markPatchPoints(boolList&);

        //- get the surface engine
        const meshSurfaceEngine& surfaceEngine() const;

        //- get the surface partitioner
        const triSurfacePartitioner& partitioner() const;

        //- get the edge classifier
        const triSurfaceClassifyEdges& edgeClassifier() const;

        //- find faces which have faces assigned to other patches
        //- as their neighbours
        void findFaceCandidates
        (
            labelLongList& faceCandidates,
            const labelList* facePatchPtr = NULL,
            const Map<label>* otherFacePatchPtr = NULL
        ) const;

        //- find patches over edges
        void findOtherFacePatchesParallel
        (
            Map<label>& otherFacePatch,
            const labelList* facePatchPtr = NULL
        ) const;

        //- find neighbour patches over edges for a boundary face
        void findNeiPatches
        (
            const label,
            const Map<label>&,
            DynList<label>&
        ) const;

        //- calculate alignment between the boundary edge and a feature edge
        //- between the specified patches
        scalar calculateAlignmentForEdge
        (
            const label beI,
            const label patch0,
            const label patch1
        ) const;

        //- calculates the deformation energy metric for the boundary edge
        //- when it would be constrained to a feature edge between
        //- the specified patches
        scalar calculateDeformationMetricForEdge
        (
            const label beI,
            const label patch0,
            const label patch1
        ) const;

        //- calculates deformation energy metric for a face
        scalar calculateDeformationMetricForFace
        (
            const label bfI,
            const DynList<label>& neiPatches,
            const label facePatch = -1
        ) const;

        //- project face centres on the nearest location at the surface mesh
        //- and assign the patch to the patch of the surface element
        void distributeBoundaryFaces();

        //- move faces into the patch with the best alignment
        bool distributeBoundaryFacesNormalAlignment();

        //- go through the boundary faces which have at least one neighour
        //- assigned to a different patch and check which of its edges
        //- are best candidates to be used as feature edges
        void findEdgeCandidates();

        //- Disallow default bitwise copy construct
        edgeExtractor(const edgeExtractor&);

        //- Disallow default bitwise assignment
        void operator=(const edgeExtractor&);

    // Private nested classes

        class faceEvaluator
        {
            // Private data
                //- const reference to the parent class
                const edgeExtractor& extractor_;

                //- calculated patches of faces over inter-processor boundaries
                Map<label> otherFacePatch_;

                //- new status of face patches after re-assigning of patches
                const labelList* newBoundaryPatchesPtr_;

                //- calculated patches of faces over inter-processor boundaries
                //- after re-assigning of patches
                Map<label>* newOtherFacePatchPtr_;

            // Private member functions
                //- calculate patches of faces over inter-processor boundaries
                void calculateNeiPatchesParallel();
                void calculateNeiPatchesParallelNewPatches();

                //- find neighbour faces over edges
                void neiFacesOverEdges(const label, DynList<label>&) const;

                //- find processors of faces over edges
                void neiFacesProcs(const label, DynList<label>&) const;

                //- calculate neighbour patches over edges of a boundary face
                void neiPatchesOverEdges
                (
                    const label bfI,
                    const labelList& fPatches,
                    const Map<label>& otherFacePatch,
                    DynList<label> &neiPatches
                ) const;

                //- evaluate new patch for a face based on the number of
                //- common edges shared between faces in other patches
                static label bestPatchTopological
                (
                    const DynList<label>& neiPatches,
                    const label currentPatch
                );

        public:

            // Constructors
                //- Construct from edgeExtractor
                faceEvaluator(const edgeExtractor&);

            // Destructor
                ~faceEvaluator();

            // Public member functions
                //- set the values for new boundary patches
                void setNewBoundaryPatches(const labelList& newBoudaryPatches);

                //- calculate neighbour patches over edges of a face
                void neiPatchesOverEdges
                (
                    const label bfI,
                    DynList<label>& neiPatches
                ) const;

                //- return the patch that the face shall be moved to such that
                //- the number of feature edges between the two patches
                //- is minimised
                label bestPatchTopological(const label bfI) const;

                //- return the best patch for a face by comparing the patch
                //- before and after modifications to face patches
                //- have been made
                //- this function checks whether the changes of face patches
                //- cause oscillations in the procedure
                label bestPatchAfterModification(const label bfI) const;
        };

        class cornerEvaluator
        {
            // Private data
                //- const reference to the parent class
                const edgeExtractor& extractor_;

                //- const reference to meshSurfacePartitioner
                const meshSurfacePartitioner& partitioner_;

                //- faces attached to points at inter-processor boundaries
                std::map<label, DynList<DynList<labelledPoint, 6> > > faceMap_;
                std::map<label, DynList<label> > facePatches_;
                std::map<label, DynList<label> > faceAtProc_;

            // Private member functions
                //- create addressing at inter-processor boundaries
                void createParallelAddressing();

                //- sort faces attached to a boundary point
                void sortedFacesAtPoint(const label, DynList<label>&) const;

        public:

            // Constructors
                //- construct from edgeExtractor and meshSurfacePartitioner
                cornerEvaluator
                (
                    const edgeExtractor&,
                    const meshSurfacePartitioner&
                );

            // Destructor
                ~cornerEvaluator();

            // Public member functions

                //- check patches of the faces attached to corners
                //- and re-assign the patches such that the mesh undergoes
                //- minimum rotation during the projection stage
                void improveCorners(labelList& newBoundaryPatches);
        };

    // Private enumerators

        enum edgeClassifier_
        {
            NONE = 0,
            SINGLECELLEDGE = 1,
            FITSATCONVEXEDGE = 2,
            FITSATCONCAVEEDGE = 4,
            ATTACHEDTOCORNER = 8,
            CANDIDATE = 16
        };

public:

    // Constructors

        //- Construct from mesh surface and octree
        edgeExtractor
        (
            polyMeshGen& mesh,
            const meshOctree& octree
        );


    // Destructor

        ~edgeExtractor();

    // Member Functions


        //- move vertices in the vicinity of discontinuities towards the
        //- the surface mesh by reducing the deviation of the mesh from the
        //- input geometry
        void moveVerticesTowardsDiscontinuities(const label nIterations = 2);

        //- check if there exist cells at concave feature edges which have more
        //- than one face at the boundary and the faces are
        //- distributed into patches at the concave edge
        bool checkConcaveEdgeCells();

        //- check and improve the distribution of mesh faces into patches
        //- in order to minimize the number of decomposed faces
        bool checkFacePatchesTopology();

        //- checks whether there exist corners which do not exist in the surface
        //- mesh, and checks whether the locations of corners in the volume mesh
        //- are near the existing counterparts in the surface mesh
        bool checkCorners();

        //- optimise distribution of mesh faces into patches
        //- in order to get better geometric quality of the mesh
        bool checkFacePatchesGeometry();

        //- find the nearest points on the surface of the volume mesh
        //- to the corners on the surface mesh
        bool findCornerCandidates();

        //- project the estimated corners and edges onto the surface mesh
        void projectDeterminedFeatureVertices();

        //- check and untangle the surface of the volume mesh
        bool untangleSurface();

        //- assemble the above functionality into a workflow
        void extractEdges();

        //- generate a surface mesh and store the created patches
        //- this is mainly intended for debugging purposes
        const triSurf* surfaceWithPatches() const;

        //- generate a surface mesh constin of facets adjacent to the requested
        //- surface point
        const triSurf* surfaceWithPatches(const label bpI) const;

        //- update mesh with selected patches
        void updateMeshPatches();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
