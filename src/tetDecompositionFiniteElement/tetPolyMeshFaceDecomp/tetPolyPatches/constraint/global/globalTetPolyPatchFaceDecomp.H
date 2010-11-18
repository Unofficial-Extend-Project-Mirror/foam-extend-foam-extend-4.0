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
    globalTetPolyPatchFaceDecomp

Description
    In parallel decompositions of FEM meshes a number of vertices may
    be shared between more that only two processors.  These vertices
    are typically located on intersections of face processor patches
    and cause all kinds of trouble, both in terms or organisation and
    communication.

    In order to simplify this issue, a processor point communication
    mechanism caters for all such points for a parallel decomposition.
    They are colected in a single global list which gets globally
    communicated (as necessary).  A processor point patch knows which
    local processor points are in the global list and where.  In order
    for this mechanism to work, the multiply-shared points need to be
    elominated from all other processor patches as well.

    Note: This patch is created automatically for parallel runs and is
    not accessible through the run-time selection mechanism as the
    user may get the construction wrong!


SourceFiles
    globalTetPolyPatchFaceDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef globalTetPolyPatchFaceDecomp_H
#define globalTetPolyPatchFaceDecomp_H

#include "tetPolyPatchFaceDecomp.H"
#include "coupledTetPolyPatchFaceDecomp.H"
#include "edgeList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class globalTetPolyPatchFaceDecomp Declaration
\*---------------------------------------------------------------------------*/

class globalTetPolyPatchFaceDecomp
:
    public coupledTetPolyPatchFaceDecomp
{
    // Private data

        //- Global point size - number of all shared points in decomposition
        label globalPointSize_;

        //- List of local multiply shared points
        labelList meshPoints_;

        //- Addressing into the global point patch
        labelList sharedPointAddr_;

        //- Global edge size - number of all shared edges in decomposition
        label globalEdgeSize_;

        //- List of local multiply shared edges
        edgeList meshEdges_;

        //- Addressing into the global point patch
        labelList sharedEdgeAddr_;

        //- List of local cut edges
        edgeList meshCutEdges_;

        //- Local cut edge multiplication mask
        scalarField meshCutEdgeMask_;

        //- Index in the boundary mesh
        label boundaryIndex_;

    // Demand driven private data

        //- Local edge indices
        mutable labelList* localEdgeIndicesPtr_;


        //- Cut edges indices, i.e. the ones originating from the
        //  points on the patch but connected to points which are
        //  internal (not on the patch)
        mutable labelList* cutEdgeIndicesPtr_;

        // Edge addressing for cut edges.
        // Edges are separated into two lists, one where the patch point
        // is the owner and one for the neighbour.  The edges are ordered
        // such that the first n add a contribution to the vertex zero and
        // so on; the breaks in the list are stored in the "start"
        // mechanism, described in lduAddressing.  This will be used to
        // calculate the cross-processor contribution of the vector-matrix
        // multiply.

            //- Cut edges owner list
            mutable labelList* cutEdgeOwnerIndicesPtr_;

            //- Cut edge owner starts
            mutable labelList* cutEdgeOwnerStartPtr_;

            //- Cut edges neighbour list
            mutable labelList* cutEdgeNeighbourIndicesPtr_;

            //- Cut edge neighbour starts
            mutable labelList* cutEdgeNeighbourStartPtr_;


            //- Doubly cut edge indices
            mutable labelList* doubleCutEdgeIndicesPtr_;

            //- Doubly cut edge owner addressing
            mutable labelList* doubleCutOwnerPtr_;

            //- Doubly cut edge neighbour addressing
            mutable labelList* doubleCutNeighbourPtr_;

            //- Owner-neighbour-double cut multiplication mask
            mutable scalarField* ownNeiDoubleMaskPtr_;


    // Private Member Functions

        //- Disallow default construct as copy
        globalTetPolyPatchFaceDecomp
        (
            const globalTetPolyPatchFaceDecomp&
        );

        //- Disallow default assignment
        void operator=(const globalTetPolyPatchFaceDecomp&);

        // Construction of demand-driven data

            //- Calculate local edge indices
            void calcLocalEdgesIndices() const;

            //- Calculate cut edge indices
            void calcCutEdgeIndices() const;

            //- Calculate cut edge addressing
            void calcCutEdgeAddressing() const;

            //- Clear cut edge addressing
            void clearCutEdgeAddressing() const;


public:

    typedef coupledTetPolyPatchFaceDecomp CoupledPointPatch;
    typedef tetPolyBoundaryMeshFaceDecomp BoundaryMesh;


    //- Runtime type information
    TypeName("global");

    // Constructors

        //- Construct from components
        globalTetPolyPatchFaceDecomp
        (
            const label globalPointSize,
            const labelList& meshPoints,
            const labelList& sharedPointAddr,
            const label globalEdgeSize,
            const edgeList& meshEdges,
            const labelList& sharedEdgeAddr,
            const edgeList& meshCutEdges,
            const scalarField& meshCutEdgeMask,
            const tetPolyBoundaryMeshFaceDecomp& bm,
            const label index
        );


    // Destructor

        virtual ~globalTetPolyPatchFaceDecomp();

    // Member functions

        //- Return name
        virtual const word& name() const
        {
            // There can only be a single patch of this type - therefore
            // its name is hard-coded.  HJ, date deleted
            return type();
        }

        //- Return size
        virtual label size() const
        {
            return meshPoints().size();
        }

        //- Return number of faces
        virtual label nFaces() const
        {
            return 0;
        }

        //- Return total number of shared points
        virtual label globalPointSize() const
        {
            return globalPointSize_;
        }

        //- Return total number of shared edges
        virtual label globalEdgeSize() const
        {
            return globalEdgeSize_;
        }

        //- Return the index of this patch in the tetPolyBoundaryMeshFaceDecomp
        virtual label index() const
        {
            return boundaryIndex_;
        }

        //- Return mesh points
        virtual const labelList& meshPoints() const
        {
            return meshPoints_;
        }

        //- Return local points.  Not implemented
        virtual const pointField& localPoints() const;

        //- Return point normals.  Not implemented
        virtual const vectorField& pointNormals() const;

        //- Face decomposition into triangle.  Not implemented
        virtual triFaceList faceTriangles(const label faceID) const;

        //- Return complete list of faces.  Not implemented
        virtual faceList triFaces() const;


        //- Return addressing into the global points list
        const labelList& sharedPointAddr() const
        {
            return sharedPointAddr_;
        }

        //- Return mesh edges
        virtual const edgeList& meshEdges() const
        {
            return meshEdges_;
        }

        //- Return addressing into the global edge list
        const labelList& sharedEdgeAddr() const
        {
            return sharedEdgeAddr_;
        }

        //- Return cut mesh edges
        const edgeList& meshCutEdges() const
        {
            return meshCutEdges_;
        }

        //- Return cut mesh edges multiplication mask
        //  Warning: this one is in the original cut edge order.
        //  For patch matrix multiplication se the reordered one
        const scalarField& meshCutEdgeMask() const
        {
            return meshCutEdgeMask_;
        }

        // Access functions for demand driven data

            //- Return list of edge indices for edges local to the patch
            virtual const labelList& localEdgeIndices() const;

            //- Return list of edge indices for cut edges
            //  (i.e. connecting points within the patch ot points outside it)
            virtual const labelList& cutEdgeIndices() const;

            // Cut edge addressing

               //- Return cut edge owner edge indices
               const labelList& cutEdgeOwnerIndices() const;

               //- Return cut edge owner edge starts
               const labelList& cutEdgeOwnerStart() const;

               //- Return cut edge neighbour edge indices
               const labelList& cutEdgeNeighbourIndices() const;

               //- Return cut edge neighbour edge starts
               const labelList& cutEdgeNeighbourStart() const;

              // Doubly cut edge addressing

                  //- Return doubly cut edge indices
                  const labelList& doubleCutEdgeIndices() const;

                  //- Return doubly cut edge owner addressing
                  //  into current patch
                  const labelList& doubleCutOwner() const;

                  //- Return doubly cut edge neighbour addressing
                  //  into current patch
                  const labelList& doubleCutNeighbour() const;

              //- Return cut edge multiplication mask
              const scalarField& ownNeiDoubleMask() const;


        // Update topology
        virtual void updateMesh();

        // Dummy functionality for new lduInteface base class
        // To be removed.  HJ, 2/Nov/2007

        // Access

            //- Return faceCell addressing
            virtual const unallocLabelList& faceCells() const
            {
                return labelList::null();
            }


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const unallocLabelList& internalData
            ) const
            {
                return tmp<labelField>(NULL);
            }

            //- Initialise interface data transfer
            virtual void initTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const
            {}

            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const
            {
                return tmp<labelField>(NULL);
            }

            //- Initialise transfer of internal field adjacent to the interface
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& iF
            ) const
            {}

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& iF
            ) const
            {
                return tmp<labelField>(NULL);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
