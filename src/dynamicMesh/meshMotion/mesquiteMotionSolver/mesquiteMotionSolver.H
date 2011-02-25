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
    mesquiteMotionSolver

Description
    Thin interface to the Mesquite Mesh Improvement library.
    Also provides surface-mesh smoothing capabilities.

SourceFiles
    mesquiteMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef mesquiteMotionSolver_H
#define mesquiteMotionSolver_H

#include "motionSolver.H"
#include "Time.H"
#include "Map.H"
#include "Switch.H"
#include "edgeList.H"
#include "pointIOField.H"
#include "MeshObject.H"
#include "HashSet.H"
#include "Mesquite_all_headers.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class mesquiteMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class mesquiteMotionSolver
:
    public motionSolver,
    public MeshObject<polyMesh, mesquiteMotionSolver>
{

    // Private data

        //- The mesh to be moved
        const polyMesh& Mesh_;

        //- Mesh characteristics [2D/3D]
        Switch twoDMesh_;

        //- Flag to determine if arrays have been initialized
        bool arraysInitialized_;

        //- Number of points
        unsigned long nPoints_;

        //- Number of elements
        unsigned long nCells_;

        //- Number of auxiliary points
        labelList nAuxPoints_;

        //- Number of auxiliary elements
        labelList nAuxCells_;

        //- Switch to toggle surface smoothing
        Switch surfaceSmoothing_;

        //- Switch to toggle volume correction
        Switch volumeCorrection_;

        //- Specify tolerance for volume correction
        scalar volCorrTolerance_;

        //- Specify max volume correction iterations
        label volCorrMaxIter_;

        // Specify tolerance for the CG solver
        scalar tolerance_;

        //- Specify multiple mesh-motion sweeps
        label nSweeps_;

        //- Specify the surface-smoothing interval
        label surfInterval_;

        //- Specify a relaxation factor for surface-smoothing
        scalar relax_;

        //- Vertex coordinate array passed into Mesquite
        mutable double* vtxCoords_;

        //- Connectivity array passed into mesquite
        mutable unsigned long* cellToNode_;

        //- Flag array for vertices (fixed/free)
        mutable int* fixFlags_;

        //- Reference points
        mutable pointIOField refPoints_;

        //- Original points prior to surface-smoothing
        mutable pointField origPoints_;

        //- The quality metric
        word qMetric_;

        //- Pointers to base element metric class
        HashTable<autoPtr<Mesquite::QualityMetric> > qMetricTable_;

        //- Pointer to base objective function
        autoPtr<Mesquite::ObjectiveFunction> objFunction_;

        //- Pointer to base optimization algorithm (vertex mover)
        autoPtr<Mesquite::VertexMover> optAlgorithm_;

        //- Termination criteria
        Mesquite::TerminationCriterion tcInner_;
        Mesquite::TerminationCriterion tcOuter_;

        //- Data specific to Laplacian surface smoothing
        labelList pIDs_;
        labelList offsets_;
        List<vectorField> pNormals_;
        List<vectorField> gradEdge_;
        List<vectorField> localPts_;
        List<scalarField> edgeMarker_;

        //- Data for the auxiliary entities
        labelList procIndices_;
        scalarField pointFractions_;
        labelListList procPointLabels_;

        List<vectorField> sendSurfFields_, recvSurfFields_;
        List<Map<label> > sendSurfPointMap_, recvSurfPointMap_;

        List<Map<label> > sendPointMap_, recvPointMap_;
        List<labelList> sendCellToNode_, recvCellToNode_;
        List<vectorField> sendPointBuffer_, recvPointBuffer_;

        //- Explicit patch coupling for surface smoothing
        Map<label> patchCoupling_;

        // Fields for the CG solver
        vectorField bV_;
        vectorField xV_;
        vectorField pV_;
        vectorField rV_;
        vectorField wV_;
        vectorField bdy_;
        scalarField pointMarker_;

        scalar oldVolume_;

    // Private Member Functions

        // Sparse Matrix multiply
        void A(const vectorField& p, vectorField& w);

        // Dot-product
        scalar dot(const vectorField& f1, const vectorField& f2);

        // CG solver
        label CG
        (
            const vectorField& b,
            vectorField& p,
            vectorField& r,
            vectorField& w,
            vectorField& x
        );

        // Compute the normalization factor for the matrix
        scalar normFactor
        (
            const vectorField& x,
            const vectorField& b,
            const vectorField& w,
            vectorField& tmpField
        );

        // Component-wise sumMag
        scalar cmptSumMag(const vectorField& field);

        // Transfer buffers after divergence compute.
        void transferBuffers(vectorField &field);

        // Apply boundary conditions
        void applyBCs(vectorField &field);

        // Private member function to compute connectivity data
        void initArrays();

        // Private member function to construct parallel connectivity data
        void initParallelConnectivity();

        // Parallel send
        void parWrite(const label toID, const label& data);

        // Parallel send (for lists)
        template <class Type>
        void parWrite(const label toID, const UList<Type>& data);

        // Parallel receive
        void parRead(const label fromID, label& data);

        // Parallel receive (for lists)
        template <class Type>
        void parRead(const label fromID, UList<Type>& data);

        // Identify coupled patches
        void identifyCoupledPatches();

        // Prepare for parallel surface smoothing
        void initParallelSurfaceSmoothing();

        // Prepare mesquite connectivity for parallel runs
        void initMesquiteParallelArrays();

        // Copy auxiliary points to/from buffers
        void copyAuxiliaryPoints(bool firstCopy);

        // Private member functions to read dictionary options
        void readOptions();

        // Apply fixed-value boundary conditions, if any.
        void applyFixedValuePatches();

        // Private member function to perform Laplacian surface smoothing
        void smoothSurfaces();

        // Compute quality of a tetrahedral cell
        scalar tetQuality
        (
            const label cIndex,
            const pointField& pField,
            bool returnVolume = false
        );

        // Private member function to check for invalid
        // cells and correct if necessary.
        void correctInvalidCells();

        // Correct global volume
        void correctGlobalVolume();

        // Enforce cylindrical constraints for slip-patches
        void enforceCylindricalConstraints();

        // Prepare point-normals with updated point positions
        void preparePointNormals();

        // Utility method to check validity of cells connected to a point.
        bool checkValidity
        (
            const vector& x,
            const labelList& jList,
            scalar& beta
        );

        // Clear out addressing
        void clearOut();

        // Update on topology change
        void update(const mapPolyMesh&);

public:

    //- Runtime type information
    TypeName("mesquiteMotionSolver");


    // Constructors

        //- Construct from polyMesh
        mesquiteMotionSolver(const polyMesh& mesh);

        mesquiteMotionSolver(const polyMesh& mesh, Istream& msData);

    // Destructor

        virtual ~mesquiteMotionSolver();


    // Member Functions

        //- Return reference to the mesh to be moved
        const polyMesh& mesh() const
        {
            return Mesh_;
        }

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Move points
        virtual bool movePoints() const;

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);

        //- Update topology (using meshObjectBase)
        virtual bool updateMesh(const mapPolyMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

