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
    tetPolyMeshCellDecomp

Description
    The tetrahedral decomposition of the mesh based on the underlying polyMesh.
    Uses minimal storage - the tet decomposition is provided on a cell-by-cell
    basis.

    This class will be used as the underlying structure for FEM mesh motion
    and therefore it needs to be able to create and store tet point-edge
    addressing for the solver.

    Note on addressing:
        The problem of addressing for Finite elements is two-fold. The
        first issue is addressing in the FEM matrix. The FEM coefficients are
        associated with edges in the tet decomposition. In order for the FOAM
        solver to work correctly, the coefficients need to be ordered in the
        upper triangular order. It is convinient to keep the order of edges
        and off-diagonal coefficients the same; therefore, the order of edges
        in the tet mesh is fixed by the ordering of points. The construction
        of the list from the edge list in polyMesh will be described in
        calcTetPolyMeshAddressing.C.

        The second issue is the construction of the FEM matrix.
        A fast method (Wavre 19/Sep/2003) is currently used:

        1) take a group of elements.  The discretisation will be done
        over a group.  A goup should be such that I can allocate a
        dense square matrix for all elements in the group.  Do NOT
        initialise arrays to zero.  As only a part of the group will
        be addressed, it needs to be zeroed out only for the entries
        it will have to use.  As the number of vertices per group will
        be variable, the local dense matrix changes for every group.

        2).  Allocate group-to-global and global-to-group arrays.  The
        second one is full of -1-s for vertices that are not in the
        group.  It will be initialised to -1 and is kept outside of
        the group.  This is done by going through all the vertices in
        the group and looking up a global-to-local index. If the index
        is -1, this global index has not been seen before: allocate
        next free local index to the vertex.

        3) For every element in the local-to-global, fidn the matrix
        row in ldu addressing.  For each element in a row, using
        global-to-local, find the local off-diagonal index and zero it
        out.  These are the only ones being used!

        4) Fill in local matrix based on the local vertex indices into
        the local dense matrix.

        5) Assemble global matrix by doing the following:
            - go through all group-to-global indices
            - for each group-to-global grab the row
            - for each non-zero element of the row, look up
            global-to-group.  This gives the x,y coordinates in the
            dense matrix and and insert it.

        6) Clear the global-to-local addressing by looking up the
        local-to-global indices and zeroing out only those entries.

SourceFiles
    tetPolyMeshCellDecomp.C
    calcTetPolyMeshCellDecompGeometry.C
    calcTetPolyMeshCellDecompAddressing.C
    calcTetPolyMeshCellDecompParPointData.C
    addParallelPointPatchCellDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef tetPolyMeshCellDecomp_H
#define tetPolyMeshCellDecomp_H

#include "GeoMesh.H"
#include "polyMesh.H"
#include "lduMesh.H"
#include "tetFemSolution.H"
#include "primitiveMesh.H"
#include "tetPolyBoundaryMeshCellDecomp.H"
#include "tetCellList.H"
#include "cellShapeList.H"
#include "FieldFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
template<class Type> class SquareMatrix;

class tetPolyMeshMapperCellDecomp;
class tetPolyMeshLduAddressingCellDecomp;

/*---------------------------------------------------------------------------*\
                     Class tetPolyMeshCellDecomp Declaration
\*---------------------------------------------------------------------------*/

class tetPolyMeshCellDecomp
:
    public GeoMesh<polyMesh>,
    public lduMesh,
    public tetFemSolution
{
    // Permanent data

        // Minimum mesh data

        //- Boundary mesh
        tetPolyBoundaryMeshCellDecomp boundary_;

        //- Offset in numbering to first cell centre
        label cellOffset_;


    // Demand-driven data

        //- Number of points
        mutable label nPoints_;

        //- Number of edges
        mutable label nEdges_;

        //- Number of tetrahedra
        mutable label nTets_;

        //- LDU addressing
        mutable tetPolyMeshLduAddressingCellDecomp* lduPtr_;

        //- Max number of points per cell
        mutable label maxNPointsForCell_;


        // Parallel mesh analysis tools

            //- List of points shared between processor patches
            mutable labelList* parPointsPtr_;

            //- List of edges shared between processor patches
            mutable edgeList* parEdgesPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPolyMeshCellDecomp(const tetPolyMeshCellDecomp&);

        //- Disallow default bitwise assignment
        void operator=(const tetPolyMeshCellDecomp&);


    // Private member functions to calculate demand driven data

        //- Check for and create a parallel point patch
        void addParallelPointPatch();

        //- Calculate points shared between parallel patches
        void calcParPointData() const;
        void clearOutParPointData() const;

        //- Clear out all data
        void clearOut() const;


public:

    // Declare name of the class and its debug switch
    ClassName("tetPolyMesh");

    typedef tetPolyMeshCellDecomp Mesh;
    typedef tetPolyBoundaryMeshCellDecomp BoundaryMesh;


    // Constructors

        //- Construct from components
        explicit tetPolyMeshCellDecomp(const polyMesh& pMesh);


    // Destructor

        ~tetPolyMeshCellDecomp();


    // Member Functions

        // Access

            //- Return the object registry
            virtual const objectRegistry& thisDb() const
            {
                return GeoMesh<polyMesh>::thisDb();
            }

            //- Return number of points in decomposition
            label nPoints() const;

            //- Return number of edges in decomposition
            label nEdges() const;

            //- Return number of tetrahedra in decomposition
            label nTets() const;

            //- Return number of cells in polyhedral mesh
            label nCells() const
            {
                return mesh_.nCells();
            }

            //- Return number of tetrahedra in decomposition for a cell
            label nTetsForCell(const label cellID) const;

            //- Return number of edges in decomposition for a face
            label nEdgesForFace(const label faceID) const;

            //- Return number of edges in decomposition connected to a
            //  given point
            label nEdgesForPoint(const label pointID) const;

            //- Return cell offset
            label cellOffset() const
            {
                return cellOffset_;
            }

            //- Return list of edge labels coming out of a point
            labelList edgesForPoint(const label pointID) const;

            //- Return points
            tmp<pointField> points() const;

            //- Return complete list of cell shapes. All are tetrahedra
            cellShapeList tetCells() const;

            //- Return reference to boundary mesh
            const tetPolyBoundaryMeshCellDecomp& boundary() const
            {
                return boundary_;
            }

            //- Return ldu addressing
            const lduAddressing& lduAddr() const;

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            virtual lduInterfacePtrsList interfaces() const
            {
                return boundary().interfaces();
            }

            //- Owner
            const unallocLabelList& owner() const
            {
                return lduAddr().lowerAddr();
            }

            //- Neighbour
            const unallocLabelList& neighbour() const
            {
                return lduAddr().upperAddr();
            }

            //- Return tetrahedral decomposition for cell
            tetCellList tets(const label cellID) const;

            //- Return max number of tets in a cell
            label maxNPointsForCell() const;

            //- Fill buffer with tet decomposition addressing
            //  Used for FEM matrix assembly.
            //  localToGlobalBuffer - sized to max number of vertices per cell
            //                        in the mesh
            //  globalToLocalBuffer - sized to total number of points in
            //                         the mesh and initialised to -1
            label addressing
            (
                const label cellID,
                labelList& localToGlobalBuffer,
                labelList& globalToLocalBuffer
            ) const;

            //- Clear global to local addressing
            void clearAddressing
            (
                const label cellID,
                const label nCellPoints,
                labelList& localToGlobalBuffer,
                labelList& globalToLocalBuffer
            ) const;

            //- Fill buffer with dot-products of shape functions
            // Used for FEM matrix assembly
            void gradNiDotGradNj
            (
                const label cellID,
                SquareMatrix<scalar>& denseMatrix,
                const labelList& globalToLocalBuffer
            ) const;

            //- Fill buffer with tensor products of shape functions
            // Used for FEM matrix assembly
            void gradNiGradNj
            (
                const label cellID,
                SquareMatrix<tensor>& denseMatrix,
                const labelList& globalToLocalBuffer
            ) const;

            //- Fill buffer with the volume integral distributed into vertices
            void volIntegral
            (
                const label cellID,
                scalarField& buffer,
                const labelList& globalToLocalBuffer
            ) const;


            // Parallel mesh analysis data

                //- Return parallel info
                const globalMeshData& globalData() const
                {
                    return mesh_.globalData();
                }

                //- Shared parallel points
                const labelList& parallelPoints() const;

                //- Shared parallel edges
                const edgeList& parallelEdges() const;


        // Edit

            //- Update mesh topology using the morph engine
            void updateMesh(const tetPolyMeshMapperCellDecomp& mapper);


    // Member Operators

        bool operator!=(const tetPolyMeshCellDecomp&) const;
        bool operator==(const tetPolyMeshCellDecomp&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
