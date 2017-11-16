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
    Foam::immersedBoundaryFvPatch

Description
    Immersed boundary FV patch

Notes
    Rewrite of the immersed boundary method for increased efficiency
    and robustness.  Main changes
    - using octree for fast search of eligible cells using indexedOctree
    - parallel communication using mapDistribute

Author
    Zeljko Tukovic
    Reorganisation and optimisation by Hrvoje Jasak

SourceFiles
    immersedBoundaryFvPatch.C
    immersedBoundaryFvPatchLeastSquaresFit.C
    immersedBoundaryFvPatchTriAddressing.C
    immersedBoundaryFvPatchSamplingWeights.C
    immersedBoundaryFvPatchTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryFvPatch_H
#define immersedBoundaryFvPatch_H

#include "fvPatch.H"
#include "immersedBoundaryPolyPatch.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dynamicLabelList.H"
#include "scalarList.H"
#include "vectorList.H"
#include "FieldFields.H"
#include "scalarMatrices.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                   Class immersedBoundaryFvPatch Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryFvPatch
:
    public fvPatch
{
    // Private data

        //- Reference to immersed boundary patch
        const immersedBoundaryPolyPatch& ibPolyPatch_;

        //- Finite volume mesh reference
        const fvMesh& mesh_;

        //- Time index for last update of mesh or moving boundary
        mutable label ibUpdateTimeIndex_;


    // Static data

        //- Fitting radius factor
        //  Cells within the radius will be used for the fitting function
        static const debug::tolerancesSwitch radiusFactor_;

        //- Fitting angle rejection factor (deg)
        //  Cells outside of the angle will not be used for the fitting function
        static const debug::tolerancesSwitch angleFactor_;

        //- Direct hit factor
        //  Fraction of cell size from IB point to IB centre which is considered
        //  a direct hit
        static const debug::tolerancesSwitch directHitFactor_;

        //- Maximum number of rows in cell-cell search
        static const debug::optimisationSwitch maxCellCellRows_;

        //- Sampling point distance factor
        //  Sampling point is located distFactor further away from the wall
        //  from the immersed boundary cell centre
        static const debug::tolerancesSwitch distFactor_;


    // Demand-driven data

        // Immersed boundary data

            //- Fluid cells indicator, marking only live fluid cells
            mutable volScalarField* gammaPtr_;

            //- Fluid cells indicator, marking live and IB cells
            mutable volScalarField* gammaExtPtr_;

            //- Fluid faces indicator, marking faces between live cells
            mutable surfaceScalarField* sGammaPtr_;

            //- List of fluid cells next to immersed boundary (IB cells)
            mutable labelList* ibCellsPtr_;

            //- List of faces for which one neighbour is an IB cell
            //  and another neighbour is a live fluid cell (IB faces)
            mutable labelList* ibFacesPtr_;

            //- List of IB cell index for each ibFace
            mutable labelList* ibFaceCellsPtr_;

            //- List of IB face flip:
            //  false if IB face points into IB cell (out of the live cell)
            //  true if IB face points into a live cell
            mutable boolList* ibFaceFlipsPtr_;

            //- List of fluid faces for which one neighbour is an IB cell
            //  and another neighbour is a dead cell (inside IB faces)
            mutable labelList* ibInsideFacesPtr_;

            //- List of internal faces in the region bounded by IB faces
            mutable labelList* ibInternalFacesPtr_;

            //- Points at the immersed boundary (IB points)
            //  nearest to the IB cell centres
            mutable vectorField* ibPointsPtr_;

            //- Normals at IB points
            mutable vectorField* ibNormalsPtr_;

            //- List of faces (triangles) which are part of IB mesh
            //  nearest to the IB cell centres
            mutable labelList* hitFacesPtr_;

            //- Direct hit flag, indicating direct hit of IB surface to
            //  IB cell centre
            mutable boolList* directHitPtr_;

            //- Sampling points for IB cells
            mutable vectorField* ibSamplingPointsPtr_;

            //- Interpolation weights for sampling weights
            mutable scalarListList* ibSamplingWeightsPtr_;

            //- Interpolation weights for sampling processor weights
            mutable scalarListList* ibSamplingProcWeightsPtr_;

            //- Interpolation addressing from IB points to tri faces
            mutable labelListList* cellsToTriAddrPtr_;

            //- Interpolation weights from IB points to tri faces
            mutable scalarListList* cellsToTriWeightsPtr_;

            //- Neighbour cells for immersed boundary cells
            //  (extended stencil)
            mutable labelListList* ibCellCellsPtr_;

            //- List of cells needed by other processors
            mutable labelListList* ibProcCellsPtr_;

            //- Processor mapDistribute
            mutable mapDistribute* ibMapPtr_;

            //- Centres of cells from other processors
            mutable vectorList* ibProcCentresPtr_;

            //- Gamma of cells from other processors (needed for sampling)
            mutable scalarList* ibProcGammaPtr_;

            //- Neighbour IB cell-proc-cell addressing
            mutable labelListList* ibCellProcCellsPtr_;

            //- Dead cells list
            mutable labelList* deadCellsPtr_;

            //- Extended dead cells list (dead cells + IB cells)
            mutable labelList* deadCellsExtPtr_;

            //- Dead faces list
            mutable labelList* deadFacesPtr_;

            //- List of live cells
            mutable labelList* liveCellsPtr_;

            //- Average IB cell sizes
            mutable scalarField* ibCellSizesPtr_;

            //- Inverse interpolation matrices for Dirichlet BC at the IB
            mutable PtrList<scalarRectangularMatrix>* invDirichletMatricesPtr_;

            //- Inverse interpolation matrices for Neumann BC at the IB
            mutable PtrList<scalarRectangularMatrix>* invNeumannMatricesPtr_;

            //- IB face area vectors
            mutable vectorField* ibSfPtr_;

            //- IB face area vector magnitudess
            mutable scalarField* ibMagSfPtr_;

            //- IB cell centre distances to IB
            mutable scalarField* ibDeltaPtr_;

            //- IB cell centre distances to IB
            mutable scalarField* ibSamplingPointDeltaPtr_;


        // Tri surface data

            //- Tri surface face area vectors
            mutable vectorField* triSfPtr_;

            //- Tri surface face labels contained in this mesh
            mutable dynamicLabelList triFacesInMesh_;


        // Cell surface search

            //- Cell octree for search in local region
            //  Used in donor cell search
            mutable indexedOctree<treeDataCell>* cellSearchPtr_;


    // Private Member Functions


        // Storage management

            //- Clear all demand-driven data
            void clearOut();

            //- Return cell set size estimate for hash sets
            label csEst() const;

            //- Return number of interpolation cells
            label nInterpolationCoeffs() const;


        // Make demand-driven data

            //- Make fluid cells indicator, marking only live fluid cells
            void makeGamma() const;

            //- Make extended fluid cells indicator, marking live and IB cells
            void makeGammaExt() const;

            //- Make fluid faces indicator
            void makeSGamma() const;

            //- Make list of cells next to immersed boundary
            void makeIbCells() const;

            //- Add corner points to IB cells list
            void addIbCornerCells() const;

            //- Make IB faces
            void makeIbFaces() const;

            //- Make tri addressing
            void makeTriAddressing() const;

            //- Make inside IB faces
            void makeIbInsideFaces() const;

            //- Make internal IB faces
            void makeIbInternalFaces() const;

            //- Make immersed boundary points and normals
            void makeIbPointsAndNormals() const;

            //- Make sampling point weights
            void makeIbSamplingWeights() const;

            //- Make extended IB cells stencils
            void makeIbCellCells() const;

            //- Make list of dead cells
            void makeDeadCells() const;

            //- Make extended list of dead cells
            void makeDeadCellsExt() const;

            //- Make list of dead faces
            void makeDeadFaces() const;

            //- Make list of live cells
            void makeLiveCells() const;

            //- Make immersed boundary cell sizes
            void makeIbCellSizes() const;

            //- Make face area vectors and magnitudes
            void makeIbSf() const;

            //- Make distance between IB cell centres
            //  and corresponding IB points
            void makeIbDelta() const;

            //- Make distance between IB cell centres
            //  and corresponding sample IB points
            void makeIbSamplingPointDelta() const;

            //- Make triangular surface face area vectors
            void makeTriSf() const;

            //- Calculate cell octree
            void calcCellSearch() const;

            //- Prepare a list of cell candidates for cell-cell addressing
            //  The list is sorted in increasing distance from pt
            void findCellCells
            (
                const vector& pt,
                const label cellID,
                labelList& cellCells
            ) const;

            //- Calc cell projection area
            scalar cellProjection(label cellID, const vector& dir) const;


        // Boundary evaluation matrices

            //- Make inverse Dirichlet interpolation matrices
            void makeInvDirichletMatrices() const;

            //- Make inverse Neumann interpolation matrices
            void makeInvNeumannMatrices() const;


protected:

    // Protected Member Functions

        //- Initialise the patches for moving points
        virtual void initMovePoints();

        //- Correct patches after moving points
        virtual void movePoints();


public:

    //- Runtime type information
    TypeName(immersedBoundaryPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        immersedBoundaryFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        );


    //- Destructor
    virtual ~immersedBoundaryFvPatch()
    {}


    // Member Functions

        // Access to immersed boundary components

            //- Return reference to immersed boundary polyPatch
            const immersedBoundaryPolyPatch& ibPolyPatch() const
            {
                return ibPolyPatch_;
            }

            //- Return immersed boundary surface mesh
            const triSurfaceMesh& ibMesh() const
            {
                return ibPolyPatch_.ibMesh();
            }

            bool internalFlow() const
            {
                return ibPolyPatch_.internalFlow();
            }

            //- Is the immersed boundary patch moving?
            bool movingIb() const
            {
                return ibPolyPatch_.movingIb();
            }


        // Immersed boundary data access

            //- Get fluid cells indicator, marking only live fluid cells
            const volScalarField& gamma() const;

            //- Get extended flud cells indicator, including live and IB cells
            const volScalarField& gammaExt() const;

            //- Get fluid faces indicator, marking faces between live cells
            const surfaceScalarField& sGamma() const;


            //- Return list of fluid cells next to immersed boundary (IB cells)
            const labelList& ibCells() const;

            //- Return list of faces for which one neighbour is an IB cell
            //  and another neighbour is a live fluid cell (IB faces)
            const labelList& ibFaces() const;

            //- Return list of IB cell index for each ibFace
            const labelList& ibFaceCells() const;

            //- Return list of IB face flip:
            //  false if IB face points into IB cell (out of the live cell)
            //  true if IB face points into a live cell
            const boolList& ibFaceFlips() const;

            //- Return list of fluid faces for which one neighbour is an
            //  IB cell and another neighbour is a dead cell (inside IB faces)
            const labelList& ibInsideFaces() const;

            //- Return list of internal faces in the region bounded by IB faces
            const labelList& ibInternalFaces() const;

            //- Return IB points
            const vectorField& ibPoints() const;

            //- Return IB normals
            const vectorField& ibNormals() const;

            //- Return list of triangles in IB mesh nearest
            //  nearest to IB cell centres
            const labelList& hitFaces() const;

            //- Return direct hit flag, indicating direct hit of IB surface to
            //  IB cell centre
            const boolList& directHit() const;

            //- Return IB sampling points
            const vectorField& ibSamplingPoints() const;

            //- Interpolation weights for sampling points
            const scalarListList& ibSamplingWeights() const;

            //- Processor interpolation weights for sampling points
            const scalarListList& ibSamplingProcWeights() const;

            //- Interpolation addressing from IB points to tri faces
            const labelListList& cellsToTriAddr() const;

            //- Interpolation weights from IB points to tri faces
            const scalarListList& cellsToTriWeights() const;

            //- Return IB cell extended stencil
            const labelListList& ibCellCells() const;

            //- Return neighbour proc cells
            //  These are cell data that other processors need from
            //  local processor
            const labelListList& ibProcCells() const;

            //- Processor mapDistribute
            const mapDistribute& ibMap() const;

            //- Return neighbour proc centres
            //  These are centres from neighbouring processors the
            //  local processor needs to receive
            const vectorList& ibProcCentres() const;

            //- Return neighbour cell addressing
            const labelListList& ibCellProcCells() const;

            //- Return neighbour proc gamma
            //  These are gamma values from neighbouring processors the
            //  local processor needs to receive
            const scalarList& ibProcGamma() const;

            //- Return dead cells
            const labelList& deadCells() const;

            //- Return extended dead cells
            const labelList& deadCellsExt() const;

            //- Return dead faces
            const labelList& deadFaces() const;

            //- Return live cells
            const labelList& liveCells() const;

            //- Return immersed boundary cell sizes
            const scalarField& ibCellSizes() const;

            //- Get inverse Dirichlet interpolation matrix
            const PtrList<scalarRectangularMatrix>&
            invDirichletMatrices() const;

            //- Get inverse Neumann interpolation matrix
            const PtrList<scalarRectangularMatrix>&
            invNeumannMatrices() const;


            //- Return IB face area vectors
            const vectorField& ibSf() const;

            //- Return IB face area vector magnitudes
            const scalarField& ibMagSf() const;

            //- Return distance to IB
            const scalarField& ibDelta() const;

            //- Return distance to IB
            const scalarField& ibSamplingPointDelta() const;

            //- Return triangular surface face area vectors
            const vectorField& triSf() const;

            //- Return triangular surface face centres
            const vectorField& triCf() const;

            //- Return labels of triangular faces which are inside the mesh
            const dynamicLabelList&  triFacesInMesh() const;


        // Donor cell search functionality

            //- Return cell octree for search in local region
            //  Only eligible cells are in the tree
            const indexedOctree<treeDataCell>& cellSearch() const;


        // Interpolation functions to and from triangular surface

            //- Collect ibPoint values: from tri face fields onto intersection
            //  points on the IB surface
            template<class Type>
            tmp<Field<Type> > toIbPoints
            (
                const Field<Type>& triValues
            ) const;

            //- Collect ibPoint values: from tri face fields onto intersection
            //  points on the IB surface with tmp
            template<class Type>
            tmp<Field<Type> > toIbPoints
            (
                const tmp<Field<Type> >& ttriValues
            ) const;

            //- triFace values: collect data from IB fields onto intersection
            //  points on the IB surface
            template<class Type>
            tmp<Field<Type> > toTriFaces
            (
                const Field<Type>& ibValues
            ) const;

            //- triFace values: collect data from IB fields onto intersection
            //  points on the IB surface with tmp
            template<class Type>
            tmp<Field<Type> > toTriFaces
            (
                const tmp<Field<Type> >& tibValues
            ) const;


        //- Interpolation functions to sampling points from mesh cell centres
        template<class Type>
        tmp<Field<Type> > toSamplingPoints
        (
            const Field<Type>& cellValues
        ) const;


        // Helper functions

            //- Renumber Field to corespond to triangular faces contained
            //  inside the mesh
            template<class Type>
            const tmp<Field<Type> > renumberField(const Field<Type>& f) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "immersedBoundaryFvPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
