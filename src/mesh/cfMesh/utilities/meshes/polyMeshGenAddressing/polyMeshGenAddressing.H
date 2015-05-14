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
    polyMeshGenAddressing

Description
    Cell-face mesh analysis engine

SourceFiles
    polyMeshGenAddressing.C
    polyMeshGenAddressingClear.C
    polyMeshGenAddressingCellCells.C
    polyMeshGenAddressingEdgeCells.C
    polyMeshGenAddressingPointCells.C
    polyMeshGenAddressingCells.C
    polyMeshGenAddressingEdgeFaces.C
    polyMeshGenAddressingPointFaces.C
    polyMeshGenAddressingCellEdges.C
    polyMeshGenAddressingPointEdges.C
    polyMeshGenAddressingPointPoints.C
    polyMeshGenAddressingEdges.C
    polyMeshGenAddressingCellCentresAndVols.C
    polyMeshGenAddressingFaceCentresAndAreas.C
    polyMeshGenAddressingUpdateGeometry.C
    polyMeshGenAddressingParallelAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenAddressing_H
#define polyMeshGenAddressing_H

#include "edgeList.H"
#include "polyMeshGenCells.H"
#include "boolList.H"
#include "HashSet.H"
#include "DynList.H"
#include "Map.H"
#include "VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class polyMeshGenAddressing Declaration
\*---------------------------------------------------------------------------*/

class polyMeshGenAddressing
{
    // Mesh data

        //- reference to the mesh
        const polyMeshGenCells& mesh_;

    // Shapes

        //- Edges
        mutable edgeList* edgesPtr_;


    // Connectivity

        //- Cell-cells
        mutable VRWGraph* ccPtr_;

        //- Edge-cells
        mutable VRWGraph* ecPtr_;

        //- Point-cells
        mutable VRWGraph* pcPtr_;

        //- Edge-faces
        mutable VRWGraph* efPtr_;

        //- Point-faces
        mutable VRWGraph* pfPtr_;

        //- Cell-edges
        mutable VRWGraph* cePtr_;

        //- Face-edges
        mutable VRWGraph* fePtr_;

        //- Point-edges
        mutable VRWGraph* pePtr_;

        //- Point-points
        mutable VRWGraph* ppPtr_;

        //- Cell-points
        mutable VRWGraph* cpPtr_;


    // Geometric data

        //- Cell centres
        mutable vectorField* cellCentresPtr_;

        //- Face centres
        mutable vectorField* faceCentresPtr_;

        //- Cell volumes
        mutable scalarField* cellVolumesPtr_;

        //- Face areas
        mutable vectorField* faceAreasPtr_;
        
    // Parallel demand driven data
        
        //- global point labels
        mutable labelLongList* globalPointLabelPtr_;
        
        //- global face labels
        mutable labelLongList* globalFaceLabelPtr_;
        
        //- global cell labels
        mutable labelLongList* globalCellLabelPtr_;
        
        //- global edge labels
        mutable labelLongList* globalEdgeLabelPtr_;
        
    // Parallel addressing
    
        //- processors containing a vertex
        mutable VRWGraph* pProcsPtr_;
        
        //- global to local point adressing
        mutable Map<label>* globalToLocalPointAddressingPtr_;
        
        //- neighbour processors sharing a point with this processor
        mutable DynList<label>* pointNeiProcsPtr_;
        
        //- processors containing an edge
        mutable VRWGraph* eProcsPtr_;
        
        //- global to local edge addressing
        mutable Map<label>* globalToLocalEdgeAddressingPtr_;
        
        //- neighbour processors sharing an edge with the current processor
        mutable DynList<label>* edgeNeiProcsPtr_;


    // Private member functions

        //- Disallow construct as copy
        polyMeshGenAddressing(const polyMeshGenAddressing&);

        //- Disallow default bitwise assignment
        void operator=(const polyMeshGenAddressing&);


    // Private member functions to calculate demand driven data

    // Topological calculations

        //- Calculate cell-cell addressing
        void calcCellCells() const;

        //- Calculate cell-edge addressing
        void calcEdgeCells() const;

        //- Calculate point-cell addressing
        void calcPointCells() const;

        //- Calculate edge-face addressing
        void calcEdgeFaces() const;
        
        //- Calculate face-edge addressing
        void calcFaceEdges() const;

        //- Calculate point-face addressing
        void calcPointFaces() const;

        //- Calculate edge list
        void calcCellEdges() const;

        //- Calculate point-edge addressing
        void calcPointEdges() const;

        //- Calculate point-point addressing
        void calcPointPoints() const;

        //- Calculate cell-point addressing
        void calcCellPoints() const;

        //- During edge calculation, a larger set of data is assembled.
        // Create and destroy as a set, using clearOutEdges()
        void calcEdges() const;
        void clearOutEdges();

    // Parallel calculations
    
        //- calculate global point labels
        void calcGlobalPointLabels() const;
        
        //- calculate global face labels
        void calcGlobalFaceLabels() const;
        
        //- calculate global cell labels
        void calcGlobalCellLabels() const;
        
        //- calculate global edge labels
        void calcGlobalEdgeLabels() const;

    // Geometrical calculations

        //- Calculate face centres and areas
        void calcFaceCentresAndAreas() const;
        void makeFaceCentresAndAreas
        (
            const pointFieldPMG& p,
            vectorField& fCtrs,
            vectorField& fAreas
        ) const;

        //- Calculate cell centres and volumes
        void calcCellCentresAndVols() const;
        void makeCellCentresAndVols
        (
            const vectorField& fCtrs,
            const vectorField& fAreas,
            vectorField& cellCtrs,
            scalarField& cellVols
        ) const;

        //- Calculate edge vectors
        void calcEdgeVectors() const;
        
    // Disallow bitwise construct
        //- Default constructor
        polyMeshGenAddressing();

public:

        // Static data

            ClassName("polyMeshGenAddressing");

/*            //- Estimated number of cells per edge
            static const unsigned cellsPerEdge_ = 4;

            //- Estimated number of cells per point
            static const unsigned cellsPerPoint_ = 8;

            //- Estimated number of faces per cell
            static const unsigned facesPerCell_ = 6;

            //- Estimated number of faces per edge
            static const unsigned facesPerEdge_ = 4;

            //- Estimated number of faces per point
            static const unsigned facesPerPoint_ = 12;

            //- Estimated number of edges per cell
            static const unsigned edgesPerCell_ = 12;

            //- Estimated number of edges per cell
            static const unsigned edgesPerFace_ = 4;

            //- Estimated number of edges per point
            static const unsigned edgesPerPoint_ = 6;

            //- Estimated number of points per cell
            static const unsigned pointsPerCell_ = 8;

            //- Estimated number of points per face
            static const unsigned pointsPerFace_ = 4;
*/
    // Constructors

        //- Construct from components
        polyMeshGenAddressing(const polyMeshGenCells& mesh);


    // Destructor

        ~polyMeshGenAddressing();


    // Member Functions

        // Primitive mesh data

            //- Return mesh edges
            const edgeList& edges() const;

        // Return mesh connectivity

            const VRWGraph& cellCells() const;
            // faceCells given as faceOwner and faceNeighbour
            const VRWGraph& edgeCells() const;
            const VRWGraph& pointCells() const;

            // faceFaces considered unnecessary
            const VRWGraph& edgeFaces() const;
            const VRWGraph& pointFaces() const;

            const VRWGraph& cellEdges() const;
            const VRWGraph& faceEdges() const;
            // edgeEdges considered unnecessary
            const VRWGraph& pointEdges() const;
            const VRWGraph& pointPoints() const;
            const VRWGraph& cellPoints() const;

        // Parallel addressing
        
            const labelLongList& globalPointLabel() const;
            const labelLongList& globalFaceLabel() const;
            const labelLongList& globalCellLabel() const;
            const labelLongList& globalEdgeLabel() const;
            
            const VRWGraph& pointAtProcs() const;
            const Map<label>& globalToLocalPointAddressing() const;
            const DynList<label>& pointNeiProcs() const;
            const VRWGraph& edgeAtProcs() const;
            const Map<label>& globalToLocalEdgeAddressing() const;
            const DynList<label>& edgeNeiProcs() const;

        // Geometric data (raw!)

            const vectorField& cellCentres() const;
            const vectorField& faceCentres() const;
            const scalarField& cellVolumes() const;
            const vectorField& faceAreas() const;

    //  Storage management

        //- Print a list of all the currently allocated mesh data
        void printAllocated() const;

        //- Clear geometry
        void clearGeom();

        //- Clear topological data
        void clearAddressing();
        
        //- Clear parallel addressing
        void clearParallelAddressing();

        //- Clear all geometry and addressing unnecessary for CFD
        void clearOut();

        //- Clear everything primitive, geometry and addressing
        void clearAll();
        
    // Update geometry data
        void updateGeometry(const boolList& changedFace);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
