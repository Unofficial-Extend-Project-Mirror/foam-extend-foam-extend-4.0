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
    ggiPolyPatch

Description
    Generalised grid interface (GGI) patch.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Contributor
    Martin Beaudoin, Hydro-Quebec, (2008)

SourceFiles
    ggiPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef ggiPolyPatch_H
#define ggiPolyPatch_H

#include "coupledPolyPatch.H"
#include "ggiInterpolation.H"
#include "transform.H"
#include "word.H"
#include "faceZone.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ggiPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class ggiPolyPatch
:
    public coupledPolyPatch
{
    // Private data

        //- Shadow patch name
        const word shadowName_;

        //- Interpolation zone name
        const word zoneName_;

        //- Use bridging to fix overlap error in interpolation
        Switch bridgeOverlap_;

        //- Shadow patch index.  Delayed evaluation for construction
        mutable label shadowIndex_;

        //- Interpolation zone index.  Delayed evaluation for construction
        mutable label zoneIndex_;

        //- Patch-to-patch interpolation
        mutable ggiZoneInterpolation* patchToPatchPtr_;

        //- Patch-to-zone addressing
        mutable labelList* zoneAddressingPtr_;

        //- Remote zone addressing: data needed for interpolation
        mutable labelList* remoteZoneAddressingPtr_;

        //- Reconstructed patch neighbour cell centres
        mutable vectorField* reconFaceCellCentresPtr_;


        // Parallel communication optimisation, stored on master processor

            //- Is the patch localised on a single processor
            // (single processor in a parallel run)?
            //  Used for parallel optimisation
            mutable bool* localParallelPtr_;

            //- List of zone faces indices received from each processor
            mutable labelListList* receiveAddrPtr_;

            //- List of zone faces indices to send to each processor
            mutable labelListList* sendAddrPtr_;


    // Private member functions

        //- Calculate patch-to-zone addressing
        virtual void calcZoneAddressing() const;

        //- Calculate remote patch-to-zone addressing
        virtual void calcRemoteZoneAddressing() const;

        //- Calculate interpolation
        virtual void calcPatchToPatch() const;


        // Geometry

            //- Calculate reconstructed cell centres
            void calcReconFaceCellCentres() const;

            //- Force calculation of transformation tensors
            virtual void calcTransforms();


        // Parallel communication optimisation, stored on master processor

            //- Calculate local parallel switch
            void calcLocalParallel() const;

            //- Calculate send and receive addressing
            void calcSendReceive() const;


            //- Return receive addressing
            const labelListList& receiveAddr() const;

            //- Return send addressing
            const labelListList& sendAddr() const;


        // Memory management

            //- Clear geometry
             void clearGeom();

            //- Clear out
            void clearOut();


protected:

    // Protected Member functions

        //- Is the GGI active? (zone and shadow present)
        bool active() const;

        //- Initialise the calculation of the patch addressing
        virtual void initAddressing();

        //- Calculate the patch addressing
        virtual void calcAddressing();

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry();

        //- Calculate the patch geometry
        virtual void calcGeometry();

        //- Initialise the patches for moving points
        virtual void initMovePoints(const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh();

        //- Update of the patch topology
        virtual void updateMesh();


public:

    //- Runtime type information
    TypeName("ggi");


    // Constructors

        //- Construct from components
        ggiPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct from components
        ggiPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& shadowName,
            const word& zoneName,
            const bool bridgeOverlap
        );

        //- Construct from dictionary
        ggiPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh&
        );

        //- Construct as copy, resetting the boundary mesh
        ggiPolyPatch
        (
            const ggiPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct as given the original patch and resetting the
        //  face list and boundary mesh information
        ggiPolyPatch
        (
            const ggiPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new ggiPolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new ggiPolyPatch
                (
                    refCast<const ggiPolyPatch>(*this),
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }


    // Destructor

        virtual ~ggiPolyPatch();


    // Member functions

        //- Return shadow patch name
        const word& shadowName() const
        {
            return shadowName_;
        }

        //- Return name of interpolation face zone
        const word& zoneName() const
        {
            return zoneName_;
        }

        //- Return shadow patch index
        label shadowIndex() const;

        //- Return zone patch index
        label zoneIndex() const;

        //- Return shadow patch
        const ggiPolyPatch& shadow() const;

        //- Return interpolation face zone
        const faceZone& zone() const;

        //- Is this the master side?
        bool master() const
        {
            return index() < shadowIndex();
        }

        //- Is this the slave side?
        bool slave() const
        {
            return !master();
        }

        //- Use bridging to fix overlap error in interpolation
        bool bridgeOverlap() const
        {
            return bridgeOverlap_;
        }

        //- Return patch-to-zone addressing
        const labelList& zoneAddressing() const;

        //- Return remote patch-to-zone addressing
        const labelList& remoteZoneAddressing() const;

        //- Is the patch localised on a single processor
        bool localParallel() const;

        //- Return reference to patch-to-patch interpolation
        //  Used only for addressing
        const ggiZoneInterpolation& patchToPatch() const;

        //- Expand face field to full zone size, including reduction
        //  New.  HJ, 12/Jun/2011
        template<class Type>
        tmp<Field<Type> > fastExpand(const Field<Type>& pf) const;

        //- Expand face field to full zone size, including reduction
        //  Obsolete.  HJ, 12/Jun/2011
        template<class Type>
        tmp<Field<Type> > expand(const Field<Type>& pf) const;

        //- Filter zone field to patch size
        //  Obsolete.  HJ, 12/Jun/2011
        template<class Type>
        tmp<Field<Type> > filter(const Field<Type>& zf) const;


        //- Interpolate face field: given field on a the shadow side,
        //  create an interpolated field on this side
        template<class Type>
        tmp<Field<Type> > interpolate(const Field<Type>& pf) const;

        template<class Type>
        tmp<Field<Type> > interpolate(const tmp<Field<Type> >& tpf) const;

        //- Bridge interpolated face field for uncovered faces
        template<class Type>
        void bridge
        (
            const Field<Type>& bridgeField,
            Field<Type>& ff
        ) const;

        //- Return reconstructed cell centres
        const vectorField& reconFaceCellCentres() const;

        //- Initialize ordering for primitivePatch. Does not
        //  refer to *this (except for name() and type() etc.)
        virtual void initOrder(const primitivePatch&) const;

        //- Return new ordering for primitivePatch.
        //  Ordering is -faceMap: for every face
        //  index of the new face -rotation: for every new face the clockwise
        //  shift of the original face. Return false if nothing changes
        //  (faceMap is identity, rotation is 0), true otherwise.
        virtual bool order
        (
            const primitivePatch&,
            labelList& faceMap,
            labelList& rotation
        ) const;

        //- Synchronise communications of ordering for primitivePatch
        //  Used in cases when no topological change happens locally,
        //  but is happening on other processors
        virtual void syncOrder() const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ggiPolyPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
