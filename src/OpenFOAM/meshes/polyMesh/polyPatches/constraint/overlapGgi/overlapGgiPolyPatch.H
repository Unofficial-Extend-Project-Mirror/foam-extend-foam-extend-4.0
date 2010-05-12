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
    overlapGgiPolyPatch

Description
    Partial overlap generalised grid interface (GGI) patch.  Master side
    remains unchanged and the slave side will be copied as needed
    to pave the master patch surface.

    This implies that the master patch has got a larger angular pitch than
    the slave and that master and slave are aligned at one edge.
    Master and slave will specify the pitch, based on which the expansion
    of the master side will be performed and used in interpolation.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    overlapGgiPolyPatch.C
    overlapGgiPolyPatchGeometry.C

\*---------------------------------------------------------------------------*/

#ifndef overlapGgiPolyPatch_H
#define overlapGgiPolyPatch_H

#include "coupledPolyPatch.H"
#include "standAlonePatch.H"
#include "overlapGgiInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class overlapGgiPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class overlapGgiPolyPatch
:
    public coupledPolyPatch
{
    // Private data

        //- Shadow patch name
        const word shadowName_;

        //- Shadow patch index.  Delayed evaluation for construction
        mutable label shadowIndex_;


        //- Rotation parameter for the overlap

            //- Rotation axis
            const vector rotationAxis_;

            //- Wedge angle
            const scalar angle_;


        // Interpolation

            //- Pointer to expanded slave patch
            mutable standAlonePatch* expandedSlavePtr_;

            //- Patch-to-expanded-patch interpolation
            mutable overlapGgiInterpolation* patchToPatchPtr_;


        //- Reconstructed patch neighbour cell centres
        mutable vectorField* reconFaceCellCentresPtr_;


    // Private member functions

        //- Calculate expanded slave patch
       void calcExpandedSlave() const;

        //- Return reference to expanded slave patch
        const standAlonePatch& expandedSlave() const;

        //- Calculate interpolation
       void calcPatchToPatch() const;

        //- Return reference to patch-to-patch interpolation
        const overlapGgiInterpolation& patchToPatch() const;

        //- Calculate reconstructed cell centres
        void calcReconFaceCellCentres() const;

        //- Force calculation of transformation tensors
        void calcTransforms();

        //- Check definition: angles and offsets
        void checkDefinition() const;

        //- Expand slave face field to full coverage
        template<class Type>
        tmp<Field<Type> > expandSlaveData(const Field<Type>& spf) const;

        //- Clear out
        void clearOut();


protected:

    // Protected Member functions

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
    TypeName("overlapGgi");


    // Constructors

        //- Construct from components
        overlapGgiPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct from components
        overlapGgiPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& shadowName,
            const vector& axis,
            const scalar angle
        );

        //- Construct from dictionary
        overlapGgiPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh&
        );

        //- Construct as copy, resetting the boundary mesh
        overlapGgiPolyPatch
        (
            const overlapGgiPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct as given the original patch and resetting the
        //  face list and boundary mesh information
        overlapGgiPolyPatch
        (
            const overlapGgiPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new overlapGgiPolyPatch(*this, bm));
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
                new overlapGgiPolyPatch
                (
                    refCast<const overlapGgiPolyPatch>(*this),
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }


    // Destructor

        virtual ~overlapGgiPolyPatch();


    // Member functions

        //- Return shadow patch name
        const word& shadowName() const
        {
            return shadowName_;
        }

        //- Return shadow patch index
        int shadowIndex() const;

        //- Return shadow patch
        const overlapGgiPolyPatch& shadow() const;

        //- Return rotation axis
        const vector& rotationAxis() const
        {
            return rotationAxis_;
        }

        //- Return wedge angle
        const scalar& angle() const
        {
            return angle_;
        }

        //- Return number of slave copies
        label nCopies() const;

        //- Is this the master side?
        bool master() const;

        //- Is this the slave side?
        bool slave() const
        {
            return !master();
        }

        //- Interpolate face field: given field on a the shadow side,
        //  create an interpolated field on this side
        template<class Type>
        tmp<Field<Type> > interpolate(const Field<Type>& pf) const;

        template<class Type>
        tmp<Field<Type> > interpolate(const tmp<Field<Type> >& tpf) const;


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
#   include "overlapGgiPolyPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
