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
    processorTetPolyPatchCellDecomp

Description
    Processor patch boundary needs to be such that the ordering of
    points in the patch is the same on both sides. Looking at the
    creation of the faces on both sides of the processor patch they
    need to be identical on both sides with the normals pointing in
    opposite directions.  This is achieved by calling the reverseFace
    function in the decomposition.  It is therefore possible to
    re-create the ordering of patch points on the slave side by
    reversing all the patch faces of the owner.  Warning: this will
    work only with the FOAM parallel decomposition tools - enforce the
    same in all other decomposition tools.

SourceFiles
    processorTetPolyPatchCellDecomp.C
    calcProcessorTetPolyPatchCellDecompAddr.C
    calcProcessorTetPolyPatchCellDecompPointAddr.C

\*---------------------------------------------------------------------------*/

#ifndef processorTetPolyPatchCellDecomp_H
#define processorTetPolyPatchCellDecomp_H

#include "coupledFaceTetPolyPatchCellDecomp.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class processorTetPolyPatchCellDecomp Declaration
\*---------------------------------------------------------------------------*/

class processorTetPolyPatchCellDecomp
:
    public coupledFaceTetPolyPatchCellDecomp
{
    // Private data

        const processorPolyPatch& procPolyPatch_;


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
        // multiply.  HJ, 13/Nov/2001

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
        processorTetPolyPatchCellDecomp
        (
            const processorTetPolyPatchCellDecomp&
        );

        //- Disallow default assignment
        void operator=(const processorTetPolyPatchCellDecomp&);


        // Construction of demand-driven data

            //- Calculate mesh points
            virtual void calcMeshPoints() const;

            //- Calculate local edge indices
            template<template<class> class FaceList>
            labelList calcProcLocalEdgesIndices
            (
                const PrimitivePatch<face, FaceList, const pointField&>& p
            ) const;

            //- Calculate cut edge indices
            void calcCutEdgeIndices() const;

            //- Calculate cut edge addressing
            void calcCutEdgeAddressing() const;

            //- Calculate local multiplication mask
            void calcOwnNeiDoubleMask() const;

            //- Clear cut edge addressing
            void clearCutEdgeAddressing() const;


public:

    typedef coupledFaceTetPolyPatchCellDecomp CoupledPointPatch;
    typedef tetPolyBoundaryMeshCellDecomp BoundaryMesh;


    //- Runtime type information
    TypeName(processorPolyPatch::typeName_());

    // Constructors

        //- Construct from components
        processorTetPolyPatchCellDecomp
        (
            const polyPatch& patch,
            const tetPolyBoundaryMeshCellDecomp& bm
        );


    // Destructor

        virtual ~processorTetPolyPatchCellDecomp();


    // Member functions

        //- Return processor number
        int myProcNo() const
        {
            return procPolyPatch_.myProcNo();
        }

        //- Return neigbour processor number
        int neighbProcNo() const
        {
            return procPolyPatch_.neighbProcNo();
        }

        //- Is this a master patch
        bool isMaster() const
        {
            return myProcNo() < neighbProcNo();
        }

        //- Is this a slave patch
        bool isSlave() const
        {
            return !isMaster();
        }

        //- Return the underlying processorPolyPatch
        const processorPolyPatch& procPolyPatch() const
        {
            return procPolyPatch_;
        }


        // Access functions for demand driven data

            //- Return pointField of points in patch.  Not implemented.
            virtual const pointField& localPoints() const;

            //- Return point unit normals.  Not implemented.
            virtual const vectorField& pointNormals() const;

            //- Return list of edge indices for edges local to the patch
            //  (i.e. connecting points within the patch)
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

              // There exists a possibility of having an internal edge for a
              // point on the processor patch which is in fact connected to
              // another point of the same patch.  This particular nastiness
              // introduces a deformation in the solution because the edge is
              // either multiplied twice or not at all.  For this purpose, the
              // offending edges need to be separated out and multiplied
              // appropriately.  This will only happen for cell tetrahedral
              // decomposition and is generally nasty.  HJ, 4/Dec/2001

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

#ifdef NoRepository
#   include "processorTetPolyPatchCellDecomp.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
